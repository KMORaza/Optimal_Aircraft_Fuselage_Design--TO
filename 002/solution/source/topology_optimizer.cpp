#include "topology_optimizer.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>

TopologyOptimizer::TopologyOptimizer(FEMSolver& fem_solver, double volume_fraction)
    : fem_solver_(fem_solver), material_(nullptr), volume_fraction_(volume_fraction),
      min_buckling_factor_(0.0), min_feature_size_(0.0), min_thickness_(0.0), max_stress_(0.0) {
    density_.resize(fem_solver_.getGeometry().getNumElements(), 1.0);
    sensitivity_.resize(density_.size(), 0.0);
}

void TopologyOptimizer::setMaterial(const Material& material) {
    material_ = &material;
}

void TopologyOptimizer::setConstraints(double min_buckling_factor, double min_feature_size, double min_thickness, double max_stress) {
    min_buckling_factor_ = min_buckling_factor;
    min_feature_size_ = min_feature_size;
    min_thickness_ = min_thickness;
    max_stress_ = max_stress;
}

void TopologyOptimizer::computeDensityStats(double& min_density, double& max_density) const {
    min_density = *std::min_element(density_.begin(), density_.end());
    max_density = *std::max_element(density_.begin(), density_.end());
}

void TopologyOptimizer::optimize(int max_iterations, double convergence_tol, Output& output) {
    std::vector<double> old_densities = density_;
    int adjustment_count = 0;
    const int max_adjustments = 20; // Increased for more attempts

    for (int iter = 0; iter < max_iterations; ++iter) {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Solve FEM
        fem_solver_.solve(density_);

        // Update skin thickness
        double avg_density = 0.0;
        for (double rho : density_) avg_density += rho;
        avg_density /= density_.size();
        double initial_thickness = 0.004; // Match main.cpp
        double new_thickness = std::max(min_thickness_, initial_thickness + avg_density * 0.002); // Allow growth
        fem_solver_.getGeometry().updateSkinThickness(new_thickness);

        // Compute sensitivities with stress and buckling penalties
        double max_stress = fem_solver_.computeMaxStress();
        double buckling_factor = fem_solver_.computeBucklingFactor();
        double stress_penalty = max_stress > max_stress_ ? (max_stress / max_stress_) : 1.0;
        double buckling_penalty = buckling_factor < min_buckling_factor_ ? 3.0 * (min_buckling_factor_ / buckling_factor) : 1.0; // Stronger penalty
        for (size_t e = 0; e < density_.size(); ++e) {
            double p = 3.0; // Increased SIMP penalization
            sensitivity_[e] = -p * std::pow(density_[e], p - 1) * fem_solver_.computeCompliance() * stress_penalty * buckling_penalty;
        }

        // Update density
        updateDensities();

        // Compute density statistics
        double min_density, max_density;
        computeDensityStats(min_density, max_density);

        // Output iteration data
        auto end_time = std::chrono::high_resolution_clock::now();
        double iter_time = std::chrono::duration<double>(end_time - start_time).count();
        double weight = computeWeight();
        double fatigue_life = fem_solver_.computeFatigueLife();
        double max_deflection = fem_solver_.computeMaxDeflection();
        double analytical_hoop = fem_solver_.computeAnalyticalHoopStress(55160.0);
        double analytical_long = fem_solver_.computeAnalyticalLongitudinalStress(55160.0);
        output.writeData({static_cast<double>(iter + 1), weight, fem_solver_.computeCompliance(), max_stress, 
                         buckling_factor, fatigue_life, max_deflection, analytical_hoop, 
                         analytical_long, min_density, max_density, iter_time});

        // Log progress every 10 iterations
        if (iter % 10 == 0) {
            std::cout << "Iteration " << iter + 1 << ": Weight = " << weight << " kg, Stress = " 
                      << max_stress / 1e6 << " MPa, Buckling = " << buckling_factor << ", Time = " << iter_time << " s\n";
        }

        // Check design validity
        std::string violation;
        if (!isValidDesign(violation)) {
            std::cout << "Invalid design at iteration " << iter + 1 << ": " << violation << "\n";
            if (++adjustment_count >= max_adjustments) {
                std::cout << "Max adjustments reached, terminating optimization.\n";
                break;
            }
            for (size_t e = 0; e < density_.size(); ++e) {
                density_[e] = std::min(1.0, std::max(0.1, density_[e] * 2.5)); // Stronger adjustment
            }
            continue;
        }

        // Check convergence
        if (checkConvergence(old_densities)) {
            std::cout << "Converged at iteration " << iter + 1 << "\n";
            break;
        }
        old_densities = density_;
        adjustment_count = 0; // Reset after valid design
    }
}

void TopologyOptimizer::updateDensities() {
    double l1 = 0.0, l2 = 1e9;
    double move = 0.4; // Increased for faster updates
    double vol = 0.0;
    double target_vol = volume_fraction_ * density_.size();

    while ((l2 - l1) / (l2 + l1) > 1e-3) {
        double lmid = 0.5 * (l1 + l2);
        vol = 0.0;
        for (size_t e = 0; e < density_.size(); ++e) {
            double rho_new = density_[e] * std::sqrt(-sensitivity_[e] / (lmid * material_->getDensity()));
            rho_new = std::max(0.1, std::min(1.0, std::max(density_[e] - move, std::min(density_[e] + move, rho_new))));
            density_[e] = rho_new;
            vol += rho_new;
        }
        if (vol > target_vol) {
            l1 = lmid;
        } else {
            l2 = lmid;
        }
    }

    // Density filter
    std::vector<double> density_filtered = density_;
    for (size_t e = 0; e < density_.size(); ++e) {
        density_[e] = std::max(0.1, std::min(1.0, density_filtered[e]));
    }
}

bool TopologyOptimizer::checkConvergence(const std::vector<double>& old_densities) {
    double max_change = 0.0;
    for (size_t e = 0; e < density_.size(); ++e) {
        max_change = std::max(max_change, std::abs(density_[e] - old_densities[e]));
    }
    return max_change < 1e-3;
}

bool TopologyOptimizer::isValidDesign(std::string& violation) const {
    double max_stress = fem_solver_.computeMaxStress();
    double buckling_factor = fem_solver_.computeBucklingFactor();
    double max_deflection = fem_solver_.computeMaxDeflection();

    if (max_stress > max_stress_) {
        violation = "Max stress " + std::to_string(max_stress / 1e6) + " MPa exceeds limit " + std::to_string(max_stress_ / 1e6) + " MPa";
        return false;
    }
    if (buckling_factor < min_buckling_factor_) {
        violation = "Buckling factor " + std::to_string(buckling_factor) + " is below " + std::to_string(min_buckling_factor_);
        return false;
    }
    if (max_deflection > fem_solver_.getMaxDeflection()) {
        violation = "Deflection " + std::to_string(max_deflection * 1000) + " mm exceeds limit " + std::to_string(fem_solver_.getMaxDeflection() * 1000) + " mm";
        return false;
    }
    violation = "";
    return true;
}

double TopologyOptimizer::computeWeight() const {
    double volume = 0.0;
    for (size_t e = 0; e < density_.size(); ++e) {
        double area = fem_solver_.getGeometry().getPanelArea();
        double thickness = fem_solver_.getGeometry().getSkinThickness();
        volume += density_[e] * area * thickness;
    }
    return volume * material_->getDensity();
}
