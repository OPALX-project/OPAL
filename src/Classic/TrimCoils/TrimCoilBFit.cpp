//
// Class TrimCoilBFit
//   General rational function fit
//   https://gitlab.psi.ch/OPAL/src/issues/157
//
// Copyright (c) 2018 - 2019, Jochem Snuverink, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "TrimCoils/TrimCoilBFit.h"

TrimCoilBFit::TrimCoilBFit(double bmax,
                           double rmin,
                           double rmax,
                           const std::vector<double>& coefnum,
                           const std::vector<double>& coefdenom,
                           const std::vector<double>& coefnumphi,
                           const std::vector<double>& coefdenomphi):
    TrimCoilFit(bmax, rmin, rmax, coefnum, coefdenom, coefnumphi, coefdenomphi)
{}

void TrimCoilBFit::doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz)
{
    // check range
    if (r < rmin_m || r > rmax_m) return;

    double btr = 0.0, dr = 0.0;
    calculateRationalFunction(RADIUS, r, btr, dr);
    double phi = 0.0, dphi = 0.0;
    calculateRationalFunction(PHI,    phi_rad, phi, dphi);

    //std::cout << "r " << r << " dr " <<  dr << std::endl;

    *bz += bmax_m * btr *  phi;
    *br += bmax_m * (dr *  phi + btr*dphi) * z;
}