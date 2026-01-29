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
#include <vector>
#include "RecursionRelation.h"
#include "DifferentialOperator.h"

namespace polynomial {

RecursionRelation::RecursionRelation(): power_m(0), highestXorder_m(10) {
    std::vector<int> poly(1, 1);
    operator_m.setPolynomial(poly, 0, 0);
}

RecursionRelation::RecursionRelation(const std::size_t &power,
                                     const std::size_t &highestXorder):
    power_m(power), highestXorder_m(highestXorder) {
    std::vector<int> poly(1, 1);
    operator_m.setPolynomial(poly, 0, 0);
    for (std::size_t i = 0; i < power_m; i++) {
        applyOperator();
    }
}

RecursionRelation::RecursionRelation(const RecursionRelation &doperator):
    operator_m(DifferentialOperator(doperator.operator_m)),
    power_m(doperator.power_m), highestXorder_m(doperator.highestXorder_m) {
}

RecursionRelation::~RecursionRelation() {
}

RecursionRelation& RecursionRelation::operator= (
                   const RecursionRelation &recursion) {
     operator_m = recursion.operator_m;
     power_m = recursion.power_m;
     highestXorder_m = recursion.highestXorder_m;
     return *this;
}

void RecursionRelation::truncate(std::size_t highestXorder) {
    highestXorder_m = highestXorder;
    operator_m.truncate(highestXorder);
}

/*  This function increases n by 1 in the differential operator used to find \n
 *  the magnetic scalar potential
 */
void RecursionRelation::applyOperator() {
/*  Differential operator has 3 terms, make three copies of current operator */
    DifferentialOperator firstTerm(operator_m);
    DifferentialOperator secondTerm(operator_m);
    DifferentialOperator thirdTerm(operator_m);
/*  Initialise polynomials
 *  p(x) = 1 - x + x^2 - ...
 *  q(x) = 1 - 2x + 3x^2 - ...
 */
    std::vector<int> p, q;
    for (std::size_t i = 0; i <= highestXorder_m; i++) {
        p.push_back(pow(-1, i));
        q.push_back(pow(-1, i) * (i + 1));
    }
/*  Differentiate first term by x, then multiply by p(x) */
    firstTerm.differentiateX();
    firstTerm.multiplyPolynomial(Polynomial(p));
/*  Differentiate second term twice by x */
    secondTerm.differentiateX();
    secondTerm.differentiateX();
/*  Differentiate third term by s twice, then multiply by q(x) */
    thirdTerm.doubleDifferentiateS();
    thirdTerm.multiplyPolynomial(Polynomial(q));
/*  Add operators to obtain final operator */
    operator_m = DifferentialOperator(firstTerm);
    operator_m.addOperator(secondTerm);
    operator_m.addOperator(thirdTerm);
}

}
