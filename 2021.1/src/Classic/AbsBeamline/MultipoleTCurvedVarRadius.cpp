/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#include <vector>
#include "gsl/gsl_sf_pow_int.h"
#include "MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/MultipoleTFunctions/CoordinateTransform.h"

using namespace endfieldmodel;

MultipoleTCurvedVarRadius::MultipoleTCurvedVarRadius(const std::string &name):
    MultipoleTBase(name),
    maxOrderX_m(10),
    varRadiusGeometry_m(1.0, 1.0, 1.0, 1.0, 1.0),
    angle_m(0.0) {
}

MultipoleTCurvedVarRadius::MultipoleTCurvedVarRadius(
                           const MultipoleTCurvedVarRadius &right):
    MultipoleTBase(right),
    maxOrderX_m(right.maxOrderX_m),
    recursion_m(right.recursion_m),
    varRadiusGeometry_m(right.varRadiusGeometry_m),
    angle_m(right.angle_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}

MultipoleTCurvedVarRadius::~MultipoleTCurvedVarRadius() {
}

ElementBase* MultipoleTCurvedVarRadius::clone() const {
    return new MultipoleTCurvedVarRadius(*this);
}

void MultipoleTCurvedVarRadius::transformCoords(Vector_t &R) {
    std::vector<double> fringeLength = getFringeLength();
    coordinatetransform::CoordinateTransform t(R[0], R[1], R[2],
                                               getLength() / 2,
                                               fringeLength[0],
                                               fringeLength[1],
                                               (getLength() / angle_m));
    std::vector<double> r = t.getTransformation();
    R[0] = r[0];
    R[1] = r[1];
    R[2] = r[2];
}

void MultipoleTCurvedVarRadius::transformBField(Vector_t &B,
                                                const Vector_t &R) {
    double length = getLength();
    std::vector<double> fringeLength = getFringeLength();
    double prefactor = (length / angle_m) *
                       (tanh((length / 2) / fringeLength[0])
                       + tanh((length / 2) / fringeLength[1]));
    double theta = fringeLength[0] * log(cosh((R[2]
                   + (length / 2)) / fringeLength[0]))
                   - fringeLength[1] * log(cosh((R[2]
                   - (length / 2)) / fringeLength[1]));
    theta /= prefactor;
    double Bx = B[0], Bs = B[2];
    B[0] = Bx * cos(theta) - Bs * sin(theta);
    B[2] = Bx * sin(theta) + Bs * cos(theta);
}

void MultipoleTCurvedVarRadius::setMaxOrder(const std::size_t &maxOrder) {
    MultipoleTBase::setMaxOrder(maxOrder);
    std::size_t N = recursion_m.size();
    while (maxOrder >= N) {
        polynomial::RecursionRelationTwo r(N, 2 * (N + maxOrderX_m + 1));
        r.resizeX(getTransMaxOrder());
        r.truncate(maxOrderX_m);
        recursion_m.push_back(r);
        N = recursion_m.size();
    }
}

double MultipoleTCurvedVarRadius::getRadius(const double &s) {
    if (getFringeDeriv(0, s) > 1.0e-12 && angle_m != 0.0) {
        return getLength() * getFringeDeriv(0, 0)
               / (getFringeDeriv(0, s) * angle_m);
    } else {
        return 1e300; // Return -1 if radius is infinite 
    }
}

double MultipoleTCurvedVarRadius::getScaleFactor(const double &x,
                                                 const double &s) {
    return (1 + x / getRadius(s));
}

double MultipoleTCurvedVarRadius::getFn(const std::size_t &n,
                                        const double &x,
                                        const double &s) {
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(0, s);
    }
    double rho = getLength() / angle_m;
    double S_0 = getFringeDeriv(0, 0);
    double y = getFringeDeriv(0, s) / (S_0 * rho);
    double func = 0.0;
    std::vector<double> fringeDerivatives;
    for (std::size_t j = 0;
         j <= recursion_m.at(n).getMaxSDerivatives();
         j++) {
        fringeDerivatives.push_back(getFringeDeriv(j, s) / (S_0 * rho));
    }
    for (std::size_t i = 0;
         i <= recursion_m.at(n).getMaxXDerivatives();
         i++) {
        double temp = 0.0;
        for (std::size_t j = 0;
             j <= recursion_m.at(n).getMaxSDerivatives();
             j++) {
            temp += recursion_m.at(n)
                    .evaluatePolynomial(x, y, i, j, fringeDerivatives)
                    * fringeDerivatives.at(j);
        }
        func += temp * getTransDeriv(i, x);
    }
    func *= gsl_sf_pow_int(-1.0, n) * S_0 * rho;
    return func;
}
