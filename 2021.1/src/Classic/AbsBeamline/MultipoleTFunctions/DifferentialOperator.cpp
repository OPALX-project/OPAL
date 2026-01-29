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
#include <vector>
#include "DifferentialOperator.h"

namespace polynomial {

DifferentialOperator::DifferentialOperator():
    xDerivatives_m(0), sDerivatives_m(0){
    std::vector<int> temp;
    temp.push_back(1);
    polynomials_m.resize(1);
    polynomials_m[0].push_back(Polynomial(temp));
}

DifferentialOperator::DifferentialOperator(const std::size_t &xDerivatives,
                                           const std::size_t &sDerivatives):
    xDerivatives_m(xDerivatives), sDerivatives_m(sDerivatives) {
    polynomials_m.resize(xDerivatives_m + 1);
    for (std::size_t i = 0; i <= xDerivatives_m; i++) {
        polynomials_m[i].resize(sDerivatives_m + 1);
    }
}

DifferentialOperator::DifferentialOperator(
                      const DifferentialOperator &doperator):
    polynomials_m(doperator.polynomials_m),
    xDerivatives_m(doperator.xDerivatives_m),
    sDerivatives_m(doperator.sDerivatives_m) {
}

DifferentialOperator::~DifferentialOperator() {
}

DifferentialOperator& DifferentialOperator::operator= (
                      const DifferentialOperator &doperator) {
    polynomials_m = doperator.polynomials_m;
    xDerivatives_m = doperator.xDerivatives_m;
    sDerivatives_m = doperator.sDerivatives_m;
    return *this;
}

void DifferentialOperator::resizeX(const std::size_t &xDerivatives) {
    std::size_t oldxDerivatives = xDerivatives_m;
    xDerivatives_m = xDerivatives;
    polynomials_m.resize(xDerivatives_m + 1);
    if (xDerivatives_m > oldxDerivatives) {
        for (std::size_t i = oldxDerivatives + 1; i <= xDerivatives_m; i++) {
            polynomials_m[i].resize(sDerivatives_m + 1);
        }
    }
}

void DifferentialOperator::resizeS(const std::size_t &sDerivatives) {
    sDerivatives_m = sDerivatives;
    for (std::size_t i = 0; i <= xDerivatives_m; i++) {
        polynomials_m[i].resize(sDerivatives_m + 1);
    }
}

void DifferentialOperator::differentiateX() {
    resizeX(xDerivatives_m + 1);
    for (std::size_t j = 0; j <= sDerivatives_m; j++) {
        std::size_t i = xDerivatives_m;
        while (i != 0) {
            i--;
            polynomials_m[i + 1][j].addPolynomial(polynomials_m[i][j]);
            polynomials_m[i][j].differentiatePolynomial();
        }
    }
}

void DifferentialOperator::doubleDifferentiateS() {
    resizeS(sDerivatives_m + 1);
    for (std::size_t i = 0; i <= xDerivatives_m; i++) {
        std::size_t j = sDerivatives_m;
        while (j != 0) {
            j--;
            polynomials_m[i][j + 1].addPolynomial(polynomials_m[i][j]);
            polynomials_m[i][j].setZero();
        }
    }
}

void DifferentialOperator::multiplyPolynomial(const Polynomial &poly) {
    for (std::size_t i = 0; i <= xDerivatives_m; i++) {
        for (std::size_t j = 0; j <= sDerivatives_m; j++) {
            polynomials_m[i][j].multiplyPolynomial(poly);
        }
    }
}

void DifferentialOperator::setPolynomial(const std::vector<int> &poly,
                                         const std::size_t &x,
                                         const std::size_t &s) {
    if (x > xDerivatives_m) {
        resizeX(x);
    }
    if (s > sDerivatives_m) {
        resizeS(s);
    }
    polynomials_m[x][s] = Polynomial(poly);
}

void DifferentialOperator::printOperator() const {
    polynomials_m[0][0].printPolynomial();
    for (std::size_t i = 0; i <= xDerivatives_m; i++) {
        for (std::size_t j = 0; j <= sDerivatives_m; j++) {
            if ((i != 0 || j != 0) &&
                (polynomials_m[i][j].getMaxXorder() != 0 ||
                 polynomials_m[i][j].getCoefficient(0) != 0)) {
                std::cout << " + (";
                polynomials_m[i][j].printPolynomial();
                std::cout << ")" << "(d/dx)^" << i << "(d/ds)^" << j;
            }
        }
    }
    std::cout << std::endl;
}

void DifferentialOperator::addOperator(const DifferentialOperator &doperator) {
    if (xDerivatives_m < doperator.xDerivatives_m) {
        resizeX(doperator.xDerivatives_m);
    }
    if (sDerivatives_m < doperator.sDerivatives_m) {
        resizeS(doperator.sDerivatives_m);
    }
    for (std::size_t i = 0; i <= doperator.xDerivatives_m; i++) {
        for (std::size_t j = 0; j <= doperator.sDerivatives_m; j++) {
            polynomials_m[i][j].addPolynomial(doperator.polynomials_m[i][j]);
        }
    }
}

bool DifferentialOperator::isPolynomialZero(const std::size_t &x,
                                            const std::size_t &s) const {
    if (x > xDerivatives_m || s > sDerivatives_m) {
        return true;
    }
    return (polynomials_m[x][s].getCoefficient(0) == 0 &&
            polynomials_m[x][s].getMaxXorder() == 0);
}

void DifferentialOperator::truncate(const std::size_t &truncateOrder) {
    for (std::size_t i = 0; i <= xDerivatives_m; i++) {
        for (std::size_t j = 0; j <= sDerivatives_m; j++) {
            polynomials_m[i][j].truncate(truncateOrder);
        }
    }
}

double DifferentialOperator::evaluatePolynomial(
                             const double &x,
                             const std::size_t &xDerivative,
                             const std::size_t &sDerivative) const {
    if (xDerivative > xDerivatives_m || sDerivative > sDerivatives_m) {
        return 0.0;
    }
    return polynomials_m[xDerivative][sDerivative].evaluatePolynomial(x);
}

}
