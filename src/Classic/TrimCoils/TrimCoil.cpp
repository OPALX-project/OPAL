//
// Class TrimCoil
//   Abstract TrimCoil class.
//
// Copyright (c) 2018 - 2019, Matthias Frey and Jochem Snuverink,
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
// and the paper
// "Matching of turn pattern measurements for cyclotrons using multiobjective optimization"
// (https://doi.org/10.1103/PhysRevAccelBeams.22.064602)
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
#include "TrimCoil.h"

#include <cmath>

#include "Physics/Units.h"
#include "Utilities/Util.h"

TrimCoil::TrimCoil(double bmax,
                   double rmin,
                   double rmax)
{
    // convert to m
    rmin_m = rmin * Units::mm2m;
    rmax_m = rmax * Units::mm2m;
    // convert to kG
    bmax_m = bmax * Units::T2kG;
}

void TrimCoil::applyField(const double r, const double z, const double phi_rad, double *br, double *bz)
{
    if (std::abs(bmax_m) < 1e-20) return;

    double phi = Util::angle_0to2pi(phi_rad);
    // check if phi is inside [phimin_m, phimax_m]
    if (phimin_m == phimax_m || Util::angleBetweenAngles(phi, phimin_m, phimax_m))
        doApplyField(r,z,phi_rad,br,bz);
}

void TrimCoil::setAzimuth(const double phimin, const double phimax)
{
    // phi convert to rad
    phimin_m = Util::angle_0to2pi(phimin * Units::deg2rad);
    phimax_m = Util::angle_0to2pi(phimax * Units::deg2rad);
}