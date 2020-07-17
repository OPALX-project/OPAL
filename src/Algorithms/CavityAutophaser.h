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

    double track(double t,
                 const double dt,
                 const double phase,
                 std::ofstream *out = NULL) const;

    const PartData &itsReference_m;
    std::shared_ptr<Component> itsCavity_m;

    Vector_t initialR_m;
    Vector_t initialP_m;

};

#endif