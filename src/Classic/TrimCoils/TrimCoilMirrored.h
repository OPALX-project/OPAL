//
// Class TrimCoilMirrored
//   Shape mirrored from TC-15 shape
//   http://accelconf.web.cern.ch/AccelConf/ipac2017/papers/thpab077.pdf
//
// Copyright (c) 2018 - 2019, Jochem Snuverink, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef TRIM_COILMIRRORED_H
#define TRIM_COILMIRRORED_H

#include "TrimCoils/TrimCoil.h"

class TrimCoilMirrored : public TrimCoil{

public:
    TrimCoilMirrored(double bmax,
                     double rmin,
                     double rmax,
                     double slope);

    virtual ~TrimCoilMirrored() { };

private:
    /// Slope in (1 / mm)
    double bslope_m;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz);

    TrimCoilMirrored() = delete;
};

#endif //TRIM_COILMIRRORED_H
