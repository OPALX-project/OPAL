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
#include <vector>
#include <iostream>
#include "RecursionRelationTwo.h"
#include "DifferentialOperatorTwo.h"
#include "TwoPolynomial.h"

namespace polynomial {

RecursionRelationTwo::RecursionRelationTwo(): power_m(0), highestXorder_m(10) {
    std::vector<int> poly2(1, 1);
    std::vector<std::vector<int>> poly1;
    poly1.push_back(poly2);
    TwoPolynomial poly(poly1);
    operator_m.setPolynomial(poly, 0, 0);
}

RecursionRelationTwo::RecursionRelationTwo(const std::size_t &power,
                                           const std::size_t &highestXorder):
    power_m(power), highestXorder_m(highestXorder) {
    std::vector<int> poly2(1, 1);
    std::vector<std::vector<int>> poly1;
    poly1.push_back(poly2);
    TwoPolynomial poly(poly1);
    operator_m.setPolynomial(poly, 0, 0);
    for (std::size_t i = 0; i < power_m; i++) {
        applyOperator();
        truncate(highestXorder);
        sortTerms();
    }
}

RecursionRelationTwo::RecursionRelationTwo(
                      const RecursionRelationTwo &doperator):
    operator_m(DifferentialOperatorTwo(doperator.operator_m)),
    power_m(doperator.power_m), highestXorder_m(doperator.highestXorder_m) {
}

RecursionRelationTwo::~RecursionRelationTwo() {
}

RecursionRelationTwo& RecursionRelationTwo::operator= (
                      const RecursionRelationTwo &recursion) {
    operator_m = recursion.operator_m;
    power_m = recursion.power_m;
    highestXorder_m = recursion.highestXorder_m;
    return *this;
}

/*  This function increases n by 1 in the differential operator used to find \n
 *  the magnetic scalar potential
 */
void RecursionRelationTwo::applyOperator() {
/*  Differential operator has 3 terms, make three copies of current operator */
    DifferentialOperatorTwo firstTerm(operator_m);
    DifferentialOperatorTwo secondTerm(operator_m);
    DifferentialOperatorTwo thirdTerm(operator_m);
/*  Initialise polynomials
 *  p(x) = s - xs^2 + x^2s^3 - ...
 *  q(x) = 1 - x + x^2 - ...
 */
    std::vector<std::vector<int>> p, q;
    p.resize(highestXorder_m + 1);
    q.resize(highestXorder_m + 1);
    for (std::size_t i = 0; i <= highestXorder_m; i++) {
        p[i].resize(highestXorder_m + 2, 0);
        q[i].resize(highestXorder_m + 1, 0);
        p[i][i + 1] = pow(-1, i);
        q[i][i] = pow(-1, i);
    }
/*  Differentiate first term by x, then multiply by p(x) */
    firstTerm.differentiateX();
    firstTerm.multiplyPolynomial(TwoPolynomial(p));
/*  Differentiate second term twice by x */
    secondTerm.differentiateX();
    secondTerm.differentiateX();
/*  Differentiate third term by s, multiply by q(x)
 *  and then differentiate by s again
 */
    thirdTerm.differentiateS();
    thirdTerm.multiplyPolynomial(TwoPolynomial(q));
    thirdTerm.differentiateS();
    thirdTerm.multiplyPolynomial(TwoPolynomial(q));
/*  Add operators to obtain final operator */
    operator_m = DifferentialOperatorTwo(firstTerm);
    operator_m.addOperator(secondTerm);
    operator_m.addOperator(thirdTerm);
}

}
