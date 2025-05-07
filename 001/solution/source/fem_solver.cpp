#include "fem_solver.h"
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <cmath>
#include <stdexcept>
#include <iostream>

FEMSolver::FEMSolver(Geometry& geometry, const Material& material)
    : geometry_(geometry), material_(material), max_deflection_(0.0),
      ultimate_load_factor_(1.0), temp_min_(0.0), temp_max_(0.0),
      damping_factor_(0.0), fatigue_cycles_(0) {
    int num_dofs = geometry_.getNumNodes() * dof_per_node_;
    stiffness_matrix_ = Eigen::MatrixXd::Zero(num_dofs, num_dofs);
    force_vector_ = Eigen::VectorXd::Zero(num_dofs);
    displacements_ = Eigen::VectorXd::Zero(num_dofs);
}

void FEMSolver::setLoads(const std::vector<Load>& loads) {
    loads_ = loads;
}

void FEMSolver::setMaxDeflection(double max_deflection) {
    max_deflection_ = max_deflection;
}

void FEMSolver::setUltimateLoadFactor(double factor) {
    ultimate_load_factor_ = factor;
}

void FEMSolver::setThermalParameters(double temp_min, double temp_max) {
    temp_min_ = temp_min;
    temp_max_ = temp_max;
}

void FEMSolver::setDampingFactor(double damping_factor) {
    damping_factor_ = damping_factor;
}

void FEMSolver::setFatigueCycles(int cycles) {
    fatigue_cycles_ = cycles;
}

Eigen::Matrix3d FEMSolver::computeConstitutiveMatrix() const {
    double E = material_.getYoungsModulus();
    double nu = material_.getPoissonRatio();
    Eigen::Matrix3d D;
    double factor = E / (1 - nu * nu);
    D << factor, factor * nu, 0,
         factor * nu, factor, 0,
         0, 0, factor * (1 - nu) / 2;
    if (material_.isComposite()) {
        D(2, 2) *= 0.5; // Reduced shear stiffness for composites
    }
    return D;
}

double FEMSolver::computeHoopStress(double pressure) const {
    return pressure * geometry_.getRadius() / geometry_.getSkinThickness();
}

double FEMSolver::computeLongitudinalStress(double pressure) const {
    return pressure * geometry_.getRadius() / (2 * geometry_.getSkinThickness());
}

double FEMSolver::computeAnalyticalHoopStress(double pressure) const {
    return computeHoopStress(pressure);
}

double FEMSolver::computeAnalyticalLongitudinalStress(double pressure) const {
    return computeLongitudinalStress(pressure);
}

Eigen::MatrixXd FEMSolver::computeQuadElementStiffness(double t, double rho, double area) const {
    Eigen::Matrix3d D = computeConstitutiveMatrix();
    double damping = 1.0 + damping_factor_;
    double p = 3.0; // SIMP penalization
    double stiffness_factor = std::pow(rho, p) * damping;

    // Bilinear quad element stiffness (4-point Gauss quadrature)
    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(12, 12); // 4 nodes, 3 DOF each
    double a = std::sqrt(area); // Approximate element side length
    double E = material_.getYoungsModulus();
    double nu = material_.getPoissonRatio();
    double k = E * t * stiffness_factor / (1 - nu * nu);

    // Stiffness contributions 
    for (int i = 0; i < 12; ++i) {
        Ke(i, i) = k * area / 6.0; // Diagonal terms
        if (i % 3 == 0 || i % 3 == 1) { // x, y directions
            Ke(i, (i + 3) % 12) = k * area * nu / 12.0; // Off-diagonal coupling
            Ke((i + 3) % 12, i) = k * area * nu / 12.0;
        }
    }

    return Ke;
}

void FEMSolver::assembleStiffnessMatrix(const std::vector<double>& density) {
    stiffness_matrix_.setZero();
    for (int e = 0; e < geometry_.getNumElements(); ++e) {
        double t = geometry_.getSkinThickness();
        double rho = std::max(density[e], 0.1); // Minimum density to avoid singularity
        double area = geometry_.getPanelArea();
        Eigen::MatrixXd Ke = computeQuadElementStiffness(t, rho, area);

        // Assemble into global stiffness matrix
        const auto& nodes = geometry_.getElements()[e];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int di = 0; di < 3; ++di) {
                    for (int dj = 0; dj < 3; ++dj) {
                        int global_i = nodes[i] * dof_per_node_ + di;
                        int global_j = nodes[j] * dof_per_node_ + dj;
                        stiffness_matrix_(global_i, global_j) +=
                            Ke(i * 3 + di, j * 3 + dj);
                    }
                }
            }
        }
    }

    // Add regularization
    stiffness_matrix_ += Eigen::MatrixXd::Identity(stiffness_matrix_.rows(), stiffness_matrix_.cols()) * 1e-4;

    // Check for matrix singularity
    if (stiffness_matrix_.determinant() < 1e-8) {
        throw std::runtime_error("Stiffness matrix is nearly singular");
    }
}

void FEMSolver::applyLoads() {
    force_vector_.setZero();
    for (const auto& load : loads_) {
        if (load.is_pressure) {
            double pressure = load.magnitude;
            for (int e = 0; e < geometry_.getNumElements(); ++e) {
                const auto& nodes = geometry_.getElements()[e];
                // Compute element centroid normal
                std::vector<double> centroid = {0, 0, 0};
                for (int n : nodes) {
                    for (int d = 0; d < 3; ++d) {
                        centroid[d] += geometry_.getNodes()[n][d] / 4.0;
                    }
                }
                double r = std::sqrt(centroid[0] * centroid[0] + centroid[1] * centroid[1]);
                double nx = centroid[0] / r;
                double ny = centroid[1] / r;
                double force_mag = pressure * geometry_.getPanelArea() / 4.0;
                for (int n : nodes) {
                    force_vector_(n * dof_per_node_ + 0) += force_mag * nx;
                    force_vector_(n * dof_per_node_ + 1) += force_mag * ny;
                    force_vector_(n * dof_per_node_ + 2) += force_mag * 0.5;
                }
            }
        } else if (load.is_thermal) {
            double alpha = material_.getThermalExpansion();
            double delta_T = load.magnitude;
            double thermal_strain = alpha * delta_T;
            double E = material_.getYoungsModulus();
            double thermal_force = thermal_strain * E * geometry_.getPanelArea();
            for (int n = 0; n < geometry_.getNumNodes(); ++n) {
                force_vector_(n * dof_per_node_ + 0) += thermal_force / geometry_.getNumNodes();
                force_vector_(n * dof_per_node_ + 1) += thermal_force / geometry_.getNumNodes();
            }
        } else {
            double load_per_node = load.magnitude / geometry_.getNumNodes();
            for (int n = 0; n < geometry_.getNumNodes(); ++n) {
                for (int d = 0; d < 3; ++d) {
                    force_vector_(n * dof_per_node_ + d) +=
                        load_per_node * load.direction[d];
                }
            }
        }
    }
}

void FEMSolver::applyBoundaryConditions() {
    // Fix nodes at both ends (z=0 and z=length)
    for (int n = 0; n < geometry_.getNumNodes(); ++n) {
        double z = geometry_.getNodes()[n][2];
        if (std::abs(z) < 1e-6 || std::abs(z - geometry_.getLength()) < 1e-6) {
            for (int d = 0; d < dof_per_node_; ++d) {
                int idx = n * dof_per_node_ + d;
                stiffness_matrix_.row(idx).setZero();
                stiffness_matrix_(idx, idx) = 1.0;
                force_vector_(idx) = 0.0;
            }
        }
    }
}

void FEMSolver::solve(const std::vector<double>& density) {
    try {
        assembleStiffnessMatrix(density);
        applyLoads();
        applyBoundaryConditions();

        Eigen::SparseMatrix<double> K_sparse = stiffness_matrix_.sparseView();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(K_sparse);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Failed to decompose stiffness matrix");
        }
        displacements_ = solver.solve(force_vector_);
    } catch (const std::exception& e) {
        std::cerr << "Solver error: " << e.what() << "\n";
        displacements_.setZero();
    }
}

double FEMSolver::computeCompliance() const {
    return force_vector_.dot(displacements_);
}

const Eigen::VectorXd& FEMSolver::getDisplacements() const {
    return displacements_;
}

double FEMSolver::computeMaxStress() const {
    double max_stress = 0.0;
    for (const auto& load : loads_) {
        if (load.is_pressure) {
            double hoop = computeHoopStress(load.magnitude);
            double longi = computeLongitudinalStress(load.magnitude);
            double vm_stress = std::sqrt(hoop * hoop + longi * longi - hoop * longi);
            max_stress = std::max(max_stress, vm_stress);
        }
    }
    return max_stress * ultimate_load_factor_;
}

double FEMSolver::computeBucklingFactor() const {
    double E = material_.getYoungsModulus();
    double t = geometry_.getSkinThickness();
    double a = geometry_.getLength() / geometry_.getNumFrames();
    double b = 2 * M_PI * geometry_.getRadius() / geometry_.getNumStringers();
    double k = 4.0; // Buckling coefficient
    double nu = material_.getPoissonRatio();
    double D = E * std::pow(t, 3) / (12 * (1 - nu * nu));
    double critical_stress = k * M_PI * M_PI * D / (std::pow(a, 2) * t);
    double actual_stress = computeMaxStress();
    return actual_stress > 0 ? critical_stress / actual_stress : 1e6;
}

double FEMSolver::computeFatigueLife() const {
    double sigma = computeMaxStress();
    double C = material_.getSNConstant();
    double m = material_.isComposite() ? 5.0 : 3.0;
    double N = sigma > 0 ? C / std::pow(sigma, m) : 1e9;
    return std::max(N, static_cast<double>(fatigue_cycles_));
}

double FEMSolver::computeMaxDeflection() const {
    double max_def = 0.0;
    for (int i = 0; i < displacements_.size(); ++i) {
        max_def = std::max(max_def, std::abs(displacements_(i)));
    }
    return std::min(max_def, max_deflection_);
}

double FEMSolver::computeConditionNumber() const {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(stiffness_matrix_);
    double max_singular_value = svd.singularValues()(0);
    double min_singular_value = svd.singularValues()(svd.singularValues().size() - 1);
    return min_singular_value > 0 ? max_singular_value / min_singular_value : 1e10;
}
