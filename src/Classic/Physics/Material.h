#ifndef MATERIAL_H
#define MATERIAL_H

#include <map>
#include <array>
#include <string>
#include <memory>
#include <cmath>

namespace Physics {
    class Material {
    public:
        enum FitCoeffs {
            A2 = 0,
            A3,
            A4,
            A5
        };

        Material(double atomicNumber,
                 double atomicMass,
                 double massDensity,
                 double radiationLength,
                 double meanExcitationEnergy,
                 std::array<double, 4> fitCoefficients):
            atomicNumber_m(atomicNumber),
            atomicMass_m(atomicMass),
            massDensity_m(massDensity),
            radiationLength_m(radiationLength / massDensity / 100),
            meanExcitationEnergy_m(meanExcitationEnergy),
            stoppingPowerFitCoefficients_m(fitCoefficients)
        { }

        double getAtomicNumber() const;         // [1]
        double getAtomicMass() const;           // [u]
        double getMassDensity() const;          // [g cm^-3]
        double getRadiationLength() const;      // [m]
        double getMeanExcitationEnergy() const; // [eV]
        double getStoppingPowerFitCoefficients(FitCoeffs n) const;

        static std::shared_ptr<Material> getMaterial(const std::string &name);
        static std::shared_ptr<Material> addMaterial(const std::string &name,
                                                     std::shared_ptr<Material> mat_ptr);
    private:
        static
        std::map<std::string, std::shared_ptr<Material> > protoTable_sm;

        const double atomicNumber_m;
        const double atomicMass_m;
        const double massDensity_m;
        const double radiationLength_m;
        const double meanExcitationEnergy_m;
        const std::array<double,4> stoppingPowerFitCoefficients_m;
    };

    inline
    double Material::getAtomicNumber() const {
        return atomicNumber_m;
    }

    inline
    double Material::getAtomicMass() const {
        return atomicMass_m;
    }

    inline
    double Material::getMassDensity() const {
        return massDensity_m;
    }

    inline
    double Material::getRadiationLength() const {
        return radiationLength_m;
    }

    inline
    double Material::getMeanExcitationEnergy() const {
        return meanExcitationEnergy_m;
    }

    inline
    double Material::getStoppingPowerFitCoefficients(Material::FitCoeffs n) const {
        return stoppingPowerFitCoefficients_m[n];
    }

    // inline
    // double Material::getMeanExcitationEnergy() const {
    //     double Z = getAtomicNumber();
    //     if (Z < 13.0)
    //         return 12 * Z + 7;

    //     return 9.76 * Z + 58.8 * std::pow(Z, -.19);
    // }
}
#endif
