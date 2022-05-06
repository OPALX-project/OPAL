//
// Class CavityAutophaser
//
// This class determines the phase of an RF cavity for which the reference particle
// is accelerated to the highest energy.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
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
                 std::ofstream *out = nullptr) const;

    const PartData &itsReference_m;
    std::shared_ptr<Component> itsCavity_m;

    Vector_t initialR_m;
    Vector_t initialP_m;

};

#endif
