//
// Class TrimCoilFit
//   Abstract TrimCoilFit class
//   General rational function fit
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
#ifndef TRIM_COILFIT_H
#define TRIM_COILFIT_H

#include "TrimCoils/TrimCoil.h"

#include <vector>

class TrimCoilFit : public TrimCoil {

public:
    TrimCoilFit(double bmax,
                double rmin,
                double rmax,
                const std::vector<double>& coefnum,
                const std::vector<double>& coefdenom,
                const std::vector<double>& coefnumphi,
                const std::vector<double>& coefdenomphi);

    virtual ~TrimCoilFit() {};

protected:
    enum PolynomType  {NUM, DENOM, NUMPHI, DENOMPHI};
    enum FunctionType {RADIUS=0, PHI=2};

    /// calculate rational function and its first derivative
    void calculateRationalFunction(FunctionType, double value, double& quot, double& der_quot) const;
    /// calculate rational function and its first and second derivative
    void calculateRationalFunction(FunctionType, double value, double& quot, double& der_quot, double& der2_quot) const;

private:
    TrimCoilFit() = delete;

    /// rational function coefficients
    std::vector<std::vector<double>> coefs;
};

#endif //TRIM_COILFIT_H
