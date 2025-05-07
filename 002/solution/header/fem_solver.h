#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include "material.h"
#include "geometry.h"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

struct Load {
    double magnitude; // N, Pa, or °C
    std::vector<double> direction; // {x, y, z} or {0, 0, 0} for pressure/thermal
    bool is_pressure; // True for pressure loads
    bool is_thermal; // True for thermal loads

    Load(double mag, const std::vector<double>& dir, bool pressure = false, bool thermal = false)
        : magnitude(mag), direction(dir), is_pressure(pressure), is_thermal(thermal) {}
};

class FEMSolver {
public:
    FEMSolver(Geometry& geometry, const Material& material);
    void setLoads(const std::vector<Load>& loads);
    void setMaxDeflection(double max_deflection);
    void setUltimateLoadFactor(double factor);
    void setThermalParameters(double temp_min, double temp_max);
    void setDampingFactor(double damping_factor);
    void setFatigueCycles(int cycles);
    void solve(const std::vector<double>& density);
    double computeCompliance() const;
    const Eigen::VectorXd& getDisplacements() const;
    double computeMaxStress() const;
    double computeBucklingFactor() const;
    double computeFatigueLife() const;
    double computeMaxDeflection() const;
    double computeAnalyticalHoopStress(double pressure) const;
    double computeAnalyticalLongitudinalStress(double pressure) const;
    double getMaxDeflection() const { return max_deflection_; }
    Geometry& getGeometry() { return geometry_; }

private:
    void assembleStiffnessMatrix(const std::vector<double>& density);
    void applyLoads();
    void applyBoundaryConditions();
    Eigen::Matrix3d computeConstitutiveMatrix() const;
    double computeHoopStress(double pressure) const;
    double computeLongitudinalStress(double pressure) const;
    Eigen::MatrixXd computeQuadElementStiffness(double t, double rho, double area) const;

    Geometry& geometry_;
    const Material& material_;
    std::vector<Load> loads_;
    double max_deflection_;
    double ultimate_load_factor_;
    double temp_min_;
    double temp_max_;
    double damping_factor_;
    int fatigue_cycles_;
    Eigen::SparseMatrix<double> stiffness_matrix_;
    Eigen::VectorXd force_vector_;
    Eigen::VectorXd displacements_;
    int dof_per_node_ = 3; // 3D displacements
};

#endif // FEM_SOLVER_H
