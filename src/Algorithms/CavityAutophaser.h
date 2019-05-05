#ifndef CAVITYAUTOPHASER
#define CAVITYAUTOPHASER

#include "AbsBeamline/Component.h"
#include "Algorithms/PartData.h"

class CavityAutophaser {
public:
    CavityAutophaser(const PartData &ref,
                     std::shared_ptr<Component> cavity);

    ~CavityAutophaser();

    double getPhaseAtMaxEnergy(const Vector_t &R,
                               const Vector_t &P,
                               double t,
                               double dt);

private:
    double guessCavityPhase(double t);
    std::pair<double, double> optimizeCavityPhase(double initialGuess,
                                                  double t,
                                                  double dt);

    double track(Vector_t R,
                 Vector_t P,
                 double t,
                 const double dt,
                 const double phase,
                 std::ofstream *out = NULL) const;
    double getEnergyMeV(const Vector_t &P);
    double getMomentum(double kineticEnergyMeV);

    const PartData &itsReference_m;
    std::shared_ptr<Component> itsCavity_m;

    Vector_t initialR_m;
    Vector_t initialP_m;

};

inline
double CavityAutophaser::getEnergyMeV(const Vector_t &P) {
    return itsReference_m.getM() * 1e-6 * (sqrt(dot(P,P) + 1) - 1);
}

inline
double CavityAutophaser::getMomentum(double kineticEnergyMeV) {
    return sqrt(std::pow(kineticEnergyMeV / (itsReference_m.getM() * 1e-6) + 1, 2) - 1);
}

#endif