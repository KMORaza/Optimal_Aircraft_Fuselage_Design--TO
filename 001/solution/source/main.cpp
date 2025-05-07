#include "fem_solver.h"
#include "topology_optimizer.h"
#include "material.h"
#include "geometry.h"
#include "output.h"
#include <iostream>

int main() {
    // Define fuselage parameters
    double radius = 2.0; // meters
    double length = 10.0; // meters
    double skin_thickness = 0.002; // 2 mm initial
    int num_frames = 15; // Frames spaced ~700 mm
    int num_stringers = 40; // Stringers around circumference
    double pressure_diff = 55160.0; // 8 psi in Pa
    double mtow = 50000.0 * 9.81; // 50,000 kg in N
    double limit_load_factor = 2.5; // Limit load factor
    double ultimate_load_factor = 1.5 * limit_load_factor; // Ultimate load factor
    double max_deflection = 0.001 * length; // 0.1% of length
    double min_thickness = 0.0015; // 1.5 mm manufacturing constraint
    double temp_min = -50.0; // °C at cruise
    double temp_max = 50.0; // °C on ground
    double damping_factor = 0.02; // Vibration damping
    int fatigue_cycles = 50000; // Pressurization cycles

    // Material: Aluminum 2024-T3
    Material aluminum("Aluminum 2024-T3", 2.7e3, 73e9, 0.33, 345e6, 280e6, 1.2e-5, 1e7);
    // Material: CFRP 
    Material cfrp("CFRP", 1.5e3, 150e9, 0.3, 600e6, 500e6, 0.5e-6, 2e7, true);

    // Initialize geometry
    Geometry geometry(radius, length, skin_thickness, num_frames, num_stringers);

    // Define load conditions
    Load static_load(mtow, {0, 0, -1}); // Weight downward
    Load pressure_load(pressure_diff, {0, 0, 0}, true); // Internal pressure
    Load dynamic_load(mtow * limit_load_factor, {0, 1, 0}); // Lateral aerodynamic load
    Load thermal_load(temp_max - temp_min, {0, 0, 0}, false, true); // Thermal load
    std::vector<Load> loads = {static_load, pressure_load, dynamic_load, thermal_load};

    // FEM solver
    FEMSolver fem_solver(geometry, aluminum);
    fem_solver.setLoads(loads);
    fem_solver.setMaxDeflection(max_deflection);
    fem_solver.setUltimateLoadFactor(ultimate_load_factor);
    fem_solver.setThermalParameters(temp_min, temp_max);
    fem_solver.setDampingFactor(damping_factor);
    fem_solver.setFatigueCycles(fatigue_cycles);

    // Initialize topology optimizer
    TopologyOptimizer optimizer(fem_solver, 0.5); // Increased volume fraction to 50%
    optimizer.setMaterial(aluminum);
    optimizer.setConstraints(0.9, 0.01, min_thickness, aluminum.getYieldStrength());

    // Initialize output
    Output output("optimization_results.csv");
    output.writeHeader({"Iteration", "Weight (kg)", "Compliance", "Max Stress (Pa)", 
                       "Min Buckling Factor", "Min Fatigue Life (cycles)", "Max Deflection (m)", 
                       "Analytical Hoop Stress (Pa)", "Analytical Long Stress (Pa)", 
                       "Condition Number", "Min Density", "Max Density"});

    // Run optimization
    std::cout << "Starting topology optimization...\n";
    optimizer.optimize(100, 1e-3, output);
    std::cout << "Optimization complete.\n";

    // Output final results
    double final_weight = optimizer.computeWeight();
    double max_stress = fem_solver.computeMaxStress();
    double buckling_factor = fem_solver.computeBucklingFactor();
    double fatigue_life = fem_solver.computeFatigueLife();
    double computed_max_deflection = fem_solver.computeMaxDeflection();
    double analytical_hoop = fem_solver.computeAnalyticalHoopStress(pressure_diff);
    double analytical_long = fem_solver.computeAnalyticalLongitudinalStress(pressure_diff);
    double condition_number = fem_solver.computeConditionNumber();
    std::cout << "Final fuselage weight: " << final_weight << " kg\n";
    std::cout << "Maximum stress: " << max_stress / 1e6 << " MPa\n";
    std::cout << "Minimum buckling factor: " << buckling_factor << "\n";
    std::cout << "Minimum fatigue life: " << fatigue_life << " cycles\n";
    std::cout << "Maximum deflection: " << computed_max_deflection * 1000 << " mm\n";
    std::cout << "Analytical hoop stress: " << analytical_hoop / 1e6 << " MPa\n";
    std::cout << "Analytical longitudinal stress: " << analytical_long / 1e6 << " MPa\n";
    std::cout << "Stiffness matrix condition number: " << condition_number << "\n";

    return 0;
}
