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


#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_pow_int.h"
#include "MultipoleTStraight.h"

using namespace endfieldmodel;

MultipoleTStraight::MultipoleTStraight(const std::string &name):
    MultipoleTBase(name),
    straightGeometry_m(getLength()) {
}

MultipoleTStraight::MultipoleTStraight(const MultipoleTStraight &right):
    MultipoleTBase(right),
    straightGeometry_m(right.straightGeometry_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}


MultipoleTStraight::~MultipoleTStraight() {
}

ElementBase* MultipoleTStraight::clone() const {
    return new MultipoleTStraight(*this);
}

void MultipoleTStraight::transformCoords(Vector_t &/*R*/) {
    //R[2] += getBoundingBoxLength();
}

void MultipoleTStraight::setMaxOrder(const std::size_t &maxOrder) {
    MultipoleTBase::setMaxOrder(maxOrder);
}

double MultipoleTStraight::getBx(const Vector_t &R) {
    double Bx = 0.0;
    for(std::size_t n = 0; n <= getMaxOrder(); n++) {
        double f_n = 0.0;
        for(std::size_t i = 0; i <= n; i++) {
            f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i + 1, R[0]) *
                   getFringeDeriv(2 * n - 2 * i, R[2]);
        }
        f_n *= gsl_sf_pow_int(-1.0, n);
        Bx += f_n * gsl_sf_pow_int(R[1], 2 * n + 1) / gsl_sf_fact(2 * n + 1);
    }
    return Bx;
}

double MultipoleTStraight::getBs(const Vector_t &R) {
    double Bs = 0.0;
    for(std::size_t n = 0; n <= getMaxOrder(); n++) {
        double f_n = 0.0;
        for(std::size_t i = 0; i <= n; i++) {
            f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i, R[0]) *
                   getFringeDeriv(2 * n - 2 * i + 1, R[2]);
        }
        f_n *= gsl_sf_pow_int(-1.0, n);
        Bs += f_n * gsl_sf_pow_int(R[1], 2 * n + 1) / gsl_sf_fact(2 * n + 1);
    }
    return Bs;
}

double MultipoleTStraight::getFn(const std::size_t &n,
                                 const double &x, 
                                 const double &s) {
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(0, s);
    }
    double f_n = 0.0;
    for (std::size_t i = 0; i <= n; i++) {
        f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i, x) *
               getFringeDeriv(2 * n - 2 * i, s);
    }
    f_n *= gsl_sf_pow_int(-1.0, n);
    return f_n;
}
 
