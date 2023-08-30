//
// Class PartData
//   PartData represents a set of reference values for use in algorithms.
//
// Copyright (c) 200x - 2023, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Algorithms/PartData.h"

#include "Utilities/LogicalError.h"

#include <cmath>


PartData::PartData(double q, double m, double momentum) {
    charge_m = q;
    mass_m = m;
    setP(momentum);
}


PartData::PartData() {
    charge_m = 1.0;
    mass_m = 0.0;
    beta_m = 1.0;
    gamma_m = 1.0e10;
}


void PartData::setP(double p) {
    if (mass_m == 0.0) {
        throw LogicalError("PartData::setP()",
                           "Particle mass must not be zero.");
    }

    if (p == 0.0) {
        throw LogicalError("PartData::setP()",
                           "Particle momentum must not be zero.");
    }

    double e = std::sqrt(p * p + mass_m * mass_m);
    beta_m = p / e;
    gamma_m = e / mass_m;
}


void PartData::setE(double energy) {
    if (energy <= mass_m) {
        throw LogicalError("PartData::setE()", "Energy should be > mass.");
    }

    gamma_m = energy / mass_m;
    //beta = std::sqrt(energy*energy - mass*mass) / energy;
    double ginv = 1.0 / gamma_m;
    beta_m = std::sqrt((1.0 - ginv) * (1.0 + ginv));
}


void PartData::setBeta(double v) {
    if (v >= 1.0) {
        throw LogicalError("PartData::setBeta()", "Beta should be < 1.");
    }

    beta_m = v;
    gamma_m = 1.0 / std::sqrt(1.0 - beta_m * beta_m);
}


void PartData::setGamma(double v) {
    if (v <= 1.0) {
        throw LogicalError("PartData::setGamma()", "Gamma should be > 1.");
    }

    gamma_m = v;
    beta_m = std::sqrt(gamma_m * gamma_m - 1.0) / gamma_m;
}

void PartData::setMomentumTolerance(double tolerance) {
    momentumTolerance_m = tolerance;
}
