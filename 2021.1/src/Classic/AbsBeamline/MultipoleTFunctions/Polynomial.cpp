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

#include <cmath>
#include <iostream>
#include "Polynomial.h"

namespace polynomial {

Polynomial::Polynomial() {
    coefficients_m.push_back(0);
    maxXorder_m = 0;
}

Polynomial::Polynomial(const std::vector<int> &coefficients) {
    if (coefficients.size() == 0) {
        coefficients_m.push_back(0);
        return;
    }
    coefficients_m = coefficients;
    maxXorder_m = coefficients_m.size() - 1;
}

Polynomial::Polynomial(const Polynomial &poly) {
    coefficients_m = poly.coefficients_m;
    maxXorder_m = coefficients_m.size() - 1;
}

Polynomial::~Polynomial() {
}

Polynomial& Polynomial::operator= (const Polynomial& poly) {
    maxXorder_m = poly.maxXorder_m;
    coefficients_m = poly.coefficients_m;
    return *this;
}

void Polynomial::differentiatePolynomial() {
    if (maxXorder_m == 0 && coefficients_m[0] == 0) {
        return;
    }
    if (maxXorder_m == 0 && coefficients_m[0] != 0) {
        coefficients_m[0] = 0;
        return;
    }
    for (std::size_t i = 0; i < maxXorder_m; i++) {
        coefficients_m[i] = coefficients_m[i + 1] * (i + 1);
    }
    coefficients_m.pop_back();
    maxXorder_m--;
}

void Polynomial::multiplyPolynomial(const Polynomial &poly) {
    if (maxXorder_m == 0 && coefficients_m[0] == 0) {
        return;
    }
    if (poly.maxXorder_m == 0 && poly.coefficients_m[0] == 0) {
        setZero();
        return;
    }
    std::vector<int> newPoly(maxXorder_m + poly.maxXorder_m + 1, 0);
    for (std::size_t i = 0; i <= maxXorder_m; i++) {
        for (std::size_t j = 0; j <= poly.maxXorder_m; j ++) {
            newPoly[i + j] = newPoly[i + j] +
            coefficients_m[i] * poly.coefficients_m[j];
        }
    }
    coefficients_m = newPoly;
    maxXorder_m = coefficients_m.size() - 1;
}

void Polynomial::printPolynomial() const {
    std::cout << coefficients_m[0];
    for (std::size_t i = 1; i <= maxXorder_m; i++) {
        if (coefficients_m[i] >= 0) {
            std::cout << " + " << coefficients_m[i] << "x^" << i;
        } else {
            std::cout << " - " << -coefficients_m[i] << "x^" << i;
        }
    }
}

void Polynomial::addPolynomial(const Polynomial &poly) {
    if (poly.maxXorder_m > maxXorder_m) {
        maxXorder_m = poly.maxXorder_m;
        coefficients_m.resize(maxXorder_m + 1, 0);
    }
    for (std::size_t i = 0; i <= poly.maxXorder_m; i++) {
        coefficients_m[i] = coefficients_m[i] + poly.coefficients_m[i];
    }
}

void Polynomial::setCoefficient(const int &coefficient,
                                const std::size_t &order) {
    if (order > maxXorder_m) {
        coefficients_m.resize(order + 1, 0);
        maxXorder_m = coefficients_m.size() - 1;
    }
    coefficients_m[order] = coefficient;
}

double Polynomial::evaluatePolynomial(const double &x) const {
    double result = 0.0;
    std::size_t i = maxXorder_m + 1;
    while (i != 0) {
        i--;
        result = result * x + coefficients_m[i];
    }
    return result;
} 

}
