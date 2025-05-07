#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>

class Material {
public:
    Material(const std::string& name, double density, double youngs_modulus,
             double poisson_ratio, double yield_strength, double fatigue_limit,
             double thermal_expansion, double sn_constant, bool is_composite = false)
        : name_(name), density_(density), youngs_modulus_(youngs_modulus),
          poisson_ratio_(poisson_ratio), yield_strength_(yield_strength),
          fatigue_limit_(fatigue_limit), thermal_expansion_(thermal_expansion),
          sn_constant_(sn_constant), is_composite_(is_composite) {}

    double getDensity() const { return density_; }
    double getYoungsModulus() const { return youngs_modulus_; }
    double getPoissonRatio() const { return poisson_ratio_; }
    double getYieldStrength() const { return yield_strength_; }
    double getFatigueLimit() const { return fatigue_limit_; }
    double getThermalExpansion() const { return thermal_expansion_; }
    double getSNConstant() const { return sn_constant_; }
    bool isComposite() const { return is_composite_; }

private:
    std::string name_;
    double density_; // kg/m^3
    double youngs_modulus_; // Pa
    double poisson_ratio_;
    double yield_strength_; // Pa
    double fatigue_limit_; // Pa
    double thermal_expansion_; // /°C
    double sn_constant_; // S-N curve constant (stress for 1 cycle)
    bool is_composite_; // True for composites
};

#endif // MATERIAL_H