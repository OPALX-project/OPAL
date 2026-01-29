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

#ifndef DIFFERENTIAL_OPERATOR_H
#define DIFFERENTIAL_OPERATOR_H

/** ---------------------------------------------------------------------
  *
  * DifferentialOperator describes a differential operator in terms of \n
  * polynomial coefficients. \n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  *  A differential operator is a linear sum of operators in the form \n
  *  @f[p(x)\frac{\partial^n}{\partial x^n}\frac{\partial^m}{\partial s^m}@f],
  *  \n and the polynomials p(x) are stored in a n by m list.
  */

#include <vector>
#include "Polynomial.h"

namespace polynomial {

class DifferentialOperator {
public:
    /** Default constructor, initialises identity operator (constant 1) */
    DifferentialOperator();
    /** Constructor, initialises operator with zero polynomials and derivatives 
     *  up to xDerivatives in x and sDerivatives in s
     *  \param xDerivatives -> Highest derivative of x
     *  \param sDerivatives -> Highest serivative in s
     */
    DifferentialOperator(const std::size_t &xDerivatives,
                         const std::size_t &sDerivatives);
    /** Copy constructor */
    DifferentialOperator(const DifferentialOperator &doperator);
    /** Assigment operator */
    DifferentialOperator& operator= (const DifferentialOperator &doperator);
    /** Destructor, does nothing */
    ~DifferentialOperator();
    /** Set highest derivative in x to xDerivatives
     *  \param xDerivatives -> Highest derivative in x
     */
    void resizeX(const std::size_t &xDerivatives);
    /** Set highest derivative in s to sDerivatives
     *  \param sDerivatives -> Highest derivative in s
     */
    void resizeS(const std::size_t &sDerivatives);
    /** Differentiate wrt x using product rule */
    void differentiateX();
    /** Differentiate wrt s twice */
    void doubleDifferentiateS();
    /** Multiplies each term with given polynomial
     *  \param poly -> Polynomial that multiplies each term
     */
    void multiplyPolynomial(const Polynomial &poly);
    /** Assign the input polynomial \n
     *  If x or s is greater than xDerivatives_m or sDerivatives_m,
     *  they are adjusted accordingly and polynomials_m is resized
     *  \param poly -> Polynomial to be assigned to the operator
     *  \param x -> Number of x-derivatives
     *  \param s -> Number of s-derivatives
     */
    void setPolynomial(const std::vector<int> &poly,
                       const std::size_t &x,
                       const std::size_t &s);
    /** Print operator, for internal debugging */
    void printOperator() const;
    /** Add the operator to Operator, term by term
     *  \param doperator -> DifferentialOperator object to be added with
     *  this DifferentialOperator object
     */
    void addOperator(const DifferentialOperator &doperator);
    /** Returns highest derivative in x */
    std::size_t getXDerivatives() const;
    /** Returns highest derivative in s */
    std::size_t getSDerivatives() const;
    /** Check if polynomial with x x-derivaties and s s-derivatives is a
     *  zero polynomial \n
     *  If x or s is larger than xDerivatives or sDerivatives true is returned
     *  \param x -> Number of x-derivatives
     *  \param s -> Number of s-derivatives
     */
    bool isPolynomialZero(const std::size_t &x, const std::size_t &s) const;
    /** Truncate all polynomials to truncateOrder \n
     *  If truncateOrder is greater than the highest power of x then
     *  zeroes are appended to the end
     *  \param truncateOrder -> Highest order of x after truncation
     */
    void truncate(const std::size_t &truncateOrder);
    /** Evaluate polynomial with xDerivative x-derivatives and
     *  sDerivative s-derivatives \n
     *  If xDerivative or sDerivative is greater than xDerivative_m
     *  or sDerivative then zero is returned
     *  \param x -> Point at which polynomial is evaluated
     *  \param xDerivative -> Number of x-derivatives
     *  \param sDerivative -> Number of s-derivatives
     */
    double evaluatePolynomial(const double &x,
                              const std::size_t &xDerivative,
                              const std::size_t &sDerivative) const;
private:
    std::vector<std::vector<Polynomial>> polynomials_m;
    std::size_t xDerivatives_m;
    std::size_t sDerivatives_m;
};

inline
    std::size_t DifferentialOperator::getXDerivatives() const {
        return xDerivatives_m;
}
inline
    std::size_t DifferentialOperator::getSDerivatives() const {
        return sDerivatives_m;
}

}

#endif
