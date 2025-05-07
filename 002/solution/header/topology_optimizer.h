#ifndef TOPOLOGY_OPTIMIZER_H
#define TOPOLOGY_OPTIMIZER_H

#include "fem_solver.h"
#include "material.h"
#include "output.h"
#include <vector>

class TopologyOptimizer {
public:
    TopologyOptimizer(FEMSolver& fem_solver, double volume_fraction);
    void setMaterial(const Material& material);
    void setConstraints(double min_buckling_factor, double min_feature_size, double min_thickness, double max_stress);
    void optimize(int max_iterations, double convergence_tol, Output& output);
    double computeWeight() const;

private:
    void updateDensities();
    bool checkConvergence(const std::vector<double>& old_densities);
    bool isValidDesign(std::string& violation) const;
    void computeDensityStats(double& min_density, double& max_density) const;

    FEMSolver& fem_solver_;
    const Material* material_;
    double volume_fraction_;
    double min_buckling_factor_;
    double min_feature_size_;
    double min_thickness_;
    double max_stress_;
    std::vector<double> density_;
    std::vector<double> sensitivity_;
};

#endif // TOPOLOGY_OPTIMIZER_H
