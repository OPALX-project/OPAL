//
// Class TrimCoilPhaseFit
//   General rational function fit of the phase shift
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
#ifndef TRIM_COILPHASEFIT_H
#define TRIM_COILPHASEFIT_H

#include "TrimCoils/TrimCoilFit.h"

#include <vector>

class TrimCoilPhaseFit : public TrimCoilFit {

public:
    TrimCoilPhaseFit(double bmax,
                     double rmin,
                     double rmax,
                     const std::vector<double>& coefnum,
                     const std::vector<double>& coefdenom,
                     const std::vector<double>& coefnumphi,
                     const std::vector<double>& coefdenomphi);

    virtual ~TrimCoilPhaseFit() { };

private:
    TrimCoilPhaseFit() = delete;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz);
};

#endif //TRIM_COILPHASEFIT_H
