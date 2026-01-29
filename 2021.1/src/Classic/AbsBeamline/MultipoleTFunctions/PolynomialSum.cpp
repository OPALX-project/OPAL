/*
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

#include <iostream>
#include <algorithm>
#include <gsl/gsl_sf_pow_int.h>
#include "PolynomialSum.h"
#include "TwoPolynomial.h"

namespace polynomial {

PolynomialSum::PolynomialSum() {
}

PolynomialSum::PolynomialSum(const TwoPolynomial &polynomial) {
    polynomialSum_m.push_back(polynomial);
}

PolynomialSum::PolynomialSum(const PolynomialSum &polynomialSum):
    polynomialSum_m(polynomialSum.polynomialSum_m) {
}

PolynomialSum::~PolynomialSum() {
}

PolynomialSum& PolynomialSum::operator= (const PolynomialSum &sum) {
    polynomialSum_m = sum.polynomialSum_m;
    return *this;
}

void PolynomialSum::differentiateX() {
    for (std::size_t i = 0; i < polynomialSum_m.size(); i++) {
        polynomialSum_m[i].differentiateX();
    }
}

void PolynomialSum::differentiateS() {
    std::size_t pSize = polynomialSum_m.size();
    for (std::size_t i = 0; i < pSize; i++) {
        std::size_t N = polynomialSum_m[i].getNdSfactors();
        for (std::size_t j = 0; j < N; j++) {
            std::vector<std::size_t> tempdS =
                                     polynomialSum_m[i].getdSfactors();
            if (tempdS[j] != 0) {
                polynomialSum_m.push_back(polynomialSum_m[i]);
                polynomialSum_m.back().multiplyConstant(tempdS[j]);
                tempdS[j] = tempdS[j] - 1;
                if ((j + 1) == N) tempdS.push_back(0);
                tempdS[j + 1] = tempdS[j + 1] + 1;
                polynomialSum_m.back().setdSfactors(tempdS);
            }
        }
        polynomialSum_m[i].differentiateS();
    }
}

void PolynomialSum::multiplyPolynomial(const TwoPolynomial &poly) {
    for (std::size_t i = 0; i < polynomialSum_m.size(); i++) {
        polynomialSum_m[i].multiplyPolynomial(poly);
    }
}

void PolynomialSum::printPolynomial() const {
    for (std::size_t i = 0; i < polynomialSum_m.size(); i++) {
        if (!polynomialSum_m[i].isZero()) {
            std::cout << " + ";
            polynomialSum_m[i].printPolynomial();
        }
    }
}

bool PolynomialSum::isPolynomialZero(const std::size_t &p) const {
    if (p < polynomialSum_m.size()) {
        return polynomialSum_m[p].isZero();
    }
    else {
        return true;
    }
}

void PolynomialSum::truncate(const std::size_t &truncateOrder) {
    for (std::size_t i = 0; i < polynomialSum_m.size(); i++) {
        polynomialSum_m[i].truncate(truncateOrder);
    }
}

void PolynomialSum::addPolynomial(const PolynomialSum &poly) {
    std::vector<TwoPolynomial> newPoly;
    newPoly.reserve(polynomialSum_m.size() + poly.polynomialSum_m.size());
    newPoly.insert(newPoly.end(),
                   polynomialSum_m.begin(),
                   polynomialSum_m.end());
    newPoly.insert(newPoly.end(),
                   poly.polynomialSum_m.begin(),
                   poly.polynomialSum_m.end());
    polynomialSum_m = newPoly;
}

std::vector<std::size_t> PolynomialSum::getdSfactors(
                                        const std::size_t &p) const {
    if (p >= polynomialSum_m.size()) {
        std::vector<std::size_t> dummy;
        return dummy;
    }
    return polynomialSum_m[p].getdSfactors();
}

void PolynomialSum::sortTerms() {
    std::size_t N = polynomialSum_m.size();
    if (N < 2) {
        return;
    }
    std::sort(polynomialSum_m.begin(), polynomialSum_m.end());
    while (polynomialSum_m[0].isZero()) {
        polynomialSum_m.erase(polynomialSum_m.begin());
        if (polynomialSum_m.size() < 2) {
            return;
        }
    }
    std::size_t i = 1;
    while (i < polynomialSum_m.size()) {
        if (polynomialSum_m[i].isZero()) {
            polynomialSum_m.erase(polynomialSum_m.begin() + i);
            continue;
        }
        if (polynomialSum_m[i - 1] == polynomialSum_m[i]) {
            polynomialSum_m[i - 1].addPolynomial(polynomialSum_m[i]);
            polynomialSum_m.erase(polynomialSum_m.begin() + i);
        } else {
            i++;
        }
    }
}

void PolynomialSum::putSumTogether(
                    const std::vector<double> &dSvalues,
                    std::vector<std::vector<double>> &finalPolynomial) const {
    finalPolynomial.resize(1);
    finalPolynomial[0].resize(1, 0.0);
    std::size_t nx = 0, ns = 0;
    for (std::size_t p = 0; p < polynomialSum_m.size(); p++) {
        if (polynomialSum_m[p].isZero()) {
            continue;
        }
        double dSfactoreval = 1.0;
        for (std::size_t q = 0;
             q < polynomialSum_m[p].dSfactors_m.size();
             q++) {
            dSfactoreval *= gsl_sf_pow_int(dSvalues[q + 1],
                            polynomialSum_m[p].dSfactors_m[q]);
        }
        if (nx < polynomialSum_m[p].maxXorder_m) {
            finalPolynomial.resize(polynomialSum_m[p].maxXorder_m + 1);
            for (std::size_t k = nx + 1;
                 k <= polynomialSum_m[p].maxXorder_m;
                 k++) {
                finalPolynomial[nx].resize(ns + 1, 0.0);
            }
            nx = polynomialSum_m[p].maxXorder_m;
        }
        if (ns < polynomialSum_m[p].maxSorder_m) {
            ns = polynomialSum_m[p].maxSorder_m;
            for (std::size_t k = 0; k <= nx; k++) {
                finalPolynomial[k].resize(ns + 1, 0.0);
            }
        }
        for (std::size_t i = 0;
             i <= polynomialSum_m[p].maxXorder_m;
             i++) {
            for (std::size_t j = 0;
                 j <= polynomialSum_m[p].maxSorder_m;
                 j++) {
                finalPolynomial[i][j] += polynomialSum_m[p]
                                               .coefficients_m[i][j] 
                                               * dSfactoreval;
            }
        }
    }
}

double PolynomialSum::evaluatePolynomial2(
                      const double &x,
                      const double &s,
                      const std::vector<double> &dSvalues) const {
    std::vector<std::vector<double>> coefficients;
    putSumTogether(dSvalues, coefficients);
    double result = 0.0;
    std::size_t i = coefficients.size();
    while (i != 0) {
        i--;
        std::size_t j = coefficients[0].size();
        double temp = 0.0;
        while (j != 0) {
            j--;
            temp = temp * s + coefficients[i][j];
        }
        result = result * x + temp;
    }
    return result;
}

}
