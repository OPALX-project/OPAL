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
#ifndef TRIM_COIL_H
#define TRIM_COIL_H

#include "Physics/Physics.h"

class TrimCoil {

public:

    TrimCoil(double bmax, double rmin, double rmax);
    /// Apply the trim coil at position r and z to Bfields br and bz
    /// Calls virtual doApplyField
    void applyField(const double r, const double z, const double phi_rad, double *br, double *bz);
    /// Set azimuthal range
    void setAzimuth(const double phimin, const double phimax);
    virtual ~TrimCoil() { };

protected:

    /// Maximum B field (kG)
    double bmax_m;
    /// Minimum radius (m)
    double rmin_m;
    /// Maximum radius (m)
    double rmax_m;
    /// Minimal azimuth (rad)
    double phimin_m = 0.0;
    /// Maximal azimuth (rad)
    double phimax_m = Physics::two_pi;

private:

    /// virtual implementation of applyField
    virtual void doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz) = 0;
};

#endif //TRIM_COIL_H
