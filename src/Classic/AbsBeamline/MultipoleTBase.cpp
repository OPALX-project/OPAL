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

#include "MultipoleTBase.h"
#include "MultipoleT.h"
#include <gsl/gsl_sf_gamma.h>

MultipoleTBase::MultipoleTBase(MultipoleT* element) :
        element_m(element) {
}

double MultipoleTBase::getBz(const Vector_t &R) {
    double Bz = 0.0;
    std::size_t n = element_m->getMaxFOrder() + 1;
    while (n != 0) {
        n--;
        Bz = Bz * R[1] * R[1] + getFn(n, R[0], R[2])
             / gsl_sf_fact(2 * n);
    }
    return Bz;
}

double MultipoleTBase::getBx(const Vector_t &R) {
    double Bx = 0.0;
    std::size_t n = element_m->getMaxFOrder() + 1;
    while (n != 0) {
        n--;
        Bx = Bx * R[1] * R[1] + element_m->getFnDerivX(n, R[0], R[2])
             / gsl_sf_fact(2 * n + 1);
    }
    Bx *= R[1];
    return Bx;
}

double MultipoleTBase::getBs(const Vector_t &R) {
    double Bs = 0.0;
    std::size_t n = element_m->getMaxFOrder() + 1;
    while (n != 0) {
        n--;
        Bs = Bs * R[1] * R[1] + element_m->getFnDerivS(n, R[0], R[2])
             / gsl_sf_fact(2 * n + 1);
    }
    Bs *= R[1] / getScaleFactor(R[0], R[2]);
    return Bs;
}
