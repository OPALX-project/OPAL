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


#include "MultipoleTCurvedConstRadius.h"
#include "gsl/gsl_sf_pow_int.h"

using namespace endfieldmodel;

MultipoleTCurvedConstRadius::MultipoleTCurvedConstRadius(
                             const std::string &name):
    MultipoleTBase(name),
    maxOrderX_m(10),
    planarArcGeometry_m(1.0, 1.0),
    angle_m(0.0) {
}

MultipoleTCurvedConstRadius::MultipoleTCurvedConstRadius(
                             const MultipoleTCurvedConstRadius &right):
    MultipoleTBase(right),
    maxOrderX_m(right.maxOrderX_m),
    recursion_m(right.recursion_m),
    planarArcGeometry_m(right.planarArcGeometry_m),
    angle_m(right.angle_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}

MultipoleTCurvedConstRadius::~MultipoleTCurvedConstRadius() {
}

ElementBase* MultipoleTCurvedConstRadius::clone() const {
    return new MultipoleTCurvedConstRadius(*this);
}

void MultipoleTCurvedConstRadius::transformCoords(Vector_t &R) {
    double radius = getLength() / angle_m;
    double alpha = atan(R[2] / (R[0] + radius ));
    if (alpha != 0.0 && angle_m != 0.0) {
        R[0] = R[2] / sin(alpha) - radius;
        R[2] = radius * alpha;// + getBoundingBoxLength();
    } else {
        //R[2] = R[2] + getBoundingBoxLength();
    }
}

void MultipoleTCurvedConstRadius::transformBField(Vector_t &B,
                                                  const Vector_t &R) {
    double theta = R[2] * angle_m / getLength();
    double Bx = B[0];
    double Bs = B[2];
    B[0] = Bx * cos(theta) - Bs * sin(theta);
    B[2] = Bx * sin(theta) + Bs * cos(theta);
}

void MultipoleTCurvedConstRadius::setMaxOrder(const std::size_t &maxOrder) {
    MultipoleTBase::setMaxOrder(maxOrder);
    std::size_t N = recursion_m.size();
    while (maxOrder >= N) {
        polynomial::RecursionRelation r(N, 2 * (N + maxOrderX_m + 1));
        r.resizeX(getTransMaxOrder());
        r.truncate(maxOrderX_m);
        recursion_m.push_back(r);
        N = recursion_m.size();
    }
}

double MultipoleTCurvedConstRadius::getRadius(const double &s) {
    if (angle_m == 0.0) {
        return -1.0;
    } else {
        return getLength() / angle_m;
    }
}

double MultipoleTCurvedConstRadius::getScaleFactor(const double &x,
                                                   const double &s) {
    return (1 + x * angle_m / getLength());
}

double MultipoleTCurvedConstRadius::getFn(const std::size_t &n,
                                          const double &x,
                                          const double &s) {
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(0, s);
    }
    double rho = getLength() / angle_m;
    double func = 0.0;
    for (std::size_t j = 0;
         j <= recursion_m.at(n).getMaxSDerivatives();
         j++) {
        double FringeDerivj = getFringeDeriv(2 * j, s);
        for (std::size_t i = 0;
             i <= recursion_m.at(n).getMaxXDerivatives();
             i++) {
            if (recursion_m.at(n).isPolynomialZero(i, j)) {
                continue;
            }
            func += (recursion_m.at(n).evaluatePolynomial(x / rho, i, j)
                    * getTransDeriv(i, x) * FringeDerivj)
                    / gsl_sf_pow_int(rho, 2 * n - i - 2 * j);
        }
    }
    func *= gsl_sf_pow_int(-1.0, n);
    return func;
}
 
