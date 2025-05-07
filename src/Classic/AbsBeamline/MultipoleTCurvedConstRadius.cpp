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
#include "MultipoleT.h"
#include "gsl/gsl_sf_pow_int.h"

MultipoleTCurvedConstRadius::MultipoleTCurvedConstRadius(MultipoleT* element):
    MultipoleTBase(element),
    planarArcGeometry_m(1.0, 1.0) {
}

void MultipoleTCurvedConstRadius::initialise() {
    planarArcGeometry_m.setElementLength(element_m->getLength());
    planarArcGeometry_m.setCurvature(element_m->getBendAngle() / element_m->getLength());
    setMaxOrder(element_m->getMaxFOrder(), element_m->getMaxXOrder());
}

void MultipoleTCurvedConstRadius::transformCoords(Vector_t &R) {
    if(element_m->getBendAngle() != 0.0) {
        double radius = element_m->getLength() / element_m->getBendAngle();
        double alpha  = atan(R[2] / (R[0] + radius));
        if (alpha != 0.0) {
            R[0] = R[2] / sin(alpha) - radius;
            R[2] = radius * alpha;
        }
    }
    R[2] += element_m->getLength() / 2.0;  // Magnet origin at the center rather than entry
}

void MultipoleTCurvedConstRadius::transformBField(Vector_t &B, const Vector_t &R) {
    double theta = R[2] * element_m->getBendAngle() / element_m->getLength();
    double Bx = B[0];
    double Bs = B[2];
    B[0] = Bx * cos(theta) - Bs * sin(theta);
    B[2] = Bx * sin(theta) + Bs * cos(theta);
}

void MultipoleTCurvedConstRadius::setMaxOrder(size_t orderZ, size_t orderX) {
    std::size_t N = recursion_m.size();
    while (orderZ >= N) {
        polynomial::RecursionRelation r(N, 2 * (N + orderX + 1));
        r.resizeX(element_m->getTransMaxOrder());
        r.truncate(orderX);
        recursion_m.push_back(r);
        N = recursion_m.size();
    }
}

Vector_t MultipoleTCurvedConstRadius::localCartesianToOpalCartesian(const Vector_t& r) {
    Vector_t result = r;
    if(element_m->getBendAngle() != 0.0) {
        double radius = element_m->getLength() / element_m->getBendAngle();
        double theta  = element_m->getBendAngle() /2.0;
        double ds = radius * sin(theta);
        double dx = radius * (1 - cos(theta));
        result[0] = -dx;
        result[2] = ds;
    }
    return result;
}

double MultipoleTCurvedConstRadius::getScaleFactor(double x, double /*s*/) {
    return (1 + x * element_m->getBendAngle() / element_m->getLength());
}

double MultipoleTCurvedConstRadius::getFn(size_t n, double x, double s) {
    if (n == 0) {
        return element_m->getTransDeriv(0, x) * element_m->getFringeDeriv(0, s);
    }
    double rho = element_m->getLength() / element_m->getBendAngle();
    double func = 0.0;
    for (std::size_t j = 0;
         j <= recursion_m.at(n).getMaxSDerivatives();
         j++) {
        double FringeDerivj = element_m->getFringeDeriv(2 * j, s);
        for (std::size_t i = 0;
             i <= recursion_m.at(n).getMaxXDerivatives();
             i++) {
            if (recursion_m.at(n).isPolynomialZero(i, j)) {
                continue;
            }
            func += (recursion_m.at(n).evaluatePolynomial(x / rho, i, j)
                    * element_m->getTransDeriv(i, x) * FringeDerivj)
                    / gsl_sf_pow_int(rho, 2 * static_cast<int>(n) - static_cast<int>(i) -
                                              2 * static_cast<int>(j));
        }
    }
    func *= gsl_sf_pow_int(-1.0, static_cast<int>(n));
    return func;
}
 
