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

#ifndef DIFFERENTIAL_OPERATOR_TWO_H
#define DIFFERENTIAL_OPERATOR_TWO_H

/** ---------------------------------------------------------------------
  *
  * DifferentialOperatorTwo describes a differential operator in terms of \n
  * polynomial coefficients. The polynomials have two variables x and S(s). \n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  *  A differential operator is a linear sum of operators in the form \n
  *  @f[p(x, S(s))\frac{\partial^n}{\partial x^n}
  *  \frac{\partial^m}{\partial s^m}@f], \n
  *  and the polynomials p(x, S(s)) are stored in a n by m list.
  */

#include <vector>
#include "PolynomialSum.h"
#include "TwoPolynomial.h"

namespace polynomial {

class DifferentialOperatorTwo {
public:
    /** Default constructor, initialises identity operator */
    DifferentialOperatorTwo();
    /** Constructor, initialises operator with zero polynomials and
      * derivatives up to xDerivatives in x and sDerivatives in s
      * \param xDerivatives -> Highest x-derivative
      * \param sDerivatives -> Highest s-derivative
      */
    DifferentialOperatorTwo(const std::size_t &xDerivatives,
                            const std::size_t &sDerivatives);
    /** Copy constructor */
    DifferentialOperatorTwo(const DifferentialOperatorTwo &doperator);
    /** Destructor */
    ~DifferentialOperatorTwo();
    /** Assigment operator */
    DifferentialOperatorTwo& operator= (
                             const DifferentialOperatorTwo &doperator);
    /** Set highest derivative in x to xDerivatives111
     *  \param xDerivatives -> Highest x-derivative
     */
    void resizeX(const std::size_t &xDerivatives);
    /** Set highest derivative in s to xDerivatives
     *  \param sDerivatives -> Highest s-derivative
     */
    void resizeS(const std::size_t &sDerivatives);
    /** Differentiate wrt x using product rule */
    void differentiateX();
    /** Differentiate wrt s using product rule */
    void differentiateS();
    /** Multiplies each term with given polynomial
     *  \param poly -> Polynomial that is multiplied with this polynomial
     */
    void multiplyPolynomial(const TwoPolynomial &poly);
    /** Assign the polynomial poly to the operator
     *  \param poly -> Polynomial to be assigned
     *  \param x -> Number of x-derivatives
     *  \param s -> Numbe of s-derivatives
     */
    void setPolynomial(const TwoPolynomial &poly,
                       const std::size_t &x,
                       const std::size_t &s);
    /** Print operator, for internal debugging */
    void printOperator() const;
    /** Add the doperator to operator, term by term
     *  \param doperator -> Differential operator to be added 
     */
    void addOperator(const DifferentialOperatorTwo &doperator);
    /** Returns highest derivative in x */
    std::size_t getXDerivatives() const;
    /** Returns highest derivative in s */
    std::size_t getSDerivatives() const;
    /** Check if polynomial with x x-derivaties and s s-derivatives is
     *  a zero polynomial \n
     *  If x or s are out of range, true is returned
     *  \param x -> Number of x-derivatives
     *  \param s -> Number of s-derivatives
     *  \param term -> Term index
     */
    bool isPolynomialZero(const std::size_t &x,
                          const std::size_t &s,
                          const std::size_t &term) const;
    /** Truncate all polynomials to truncateOrder
     *  \param truncateOrder -> Highest order of x after truncation
     */
    void truncate(const std::size_t &truncateOrder);
    /** Evaluate polynomial specified
     *  \param x -> Point x where polynomial is evaluated
     *  \param s -> Point s where polynomial is evaluated
     *  \param xDerivative -> Number of x-derivatives
     *  \param sDerivative -> Number of s-derivatives
     *  \param term -> Term index
     */
    double evaluatePolynomial(const double &x,
                              const double &s,
                              const std::size_t &xDerivative,
                              const std::size_t &sDerivative,
                              const std::vector<double> &dSvalues) const;
    /** Returns number of terms in the sum with xDerivatives
     *  x-derivatives and sDerivatives s-derivatives \n
     *  If xDerivative or sDerivative are out of
     *  range zero is returned
     *  \param xDerivatives -> Number of x-derivatives
     *  \param sDerivatives -> Number of s-derivatives
     */
    std::size_t numberOfTerms(const std::size_t &xDerivatives,
                              const std::size_t &sDerivatives) const;
    /** Returns list of S(s)-derivatives from term p in polynomial with
     *  xDerivatives x-derivatives and sDerivatives s-derivatives \n
     *  If xDerivatives or sDerivatives are out of range
     *  an empty vector is returned
     * \param xDerivatives -> Number of x-derivatives
     * \param sDerivatives -> Number of s-derivatives
     * \param p -> Term index, starting with term 0
     */
    std::vector<std::size_t> getdSFactors(const std::size_t &xDerivatives,
                                          const std::size_t &sDerivatives,
                                          const std::size_t &p) const;
    /** Sort the terms in each sum, with fewest S(s)-derivatives first,
     *  then in increasing powers
     */
    void sortTerms();
private:
    std::vector<std::vector<PolynomialSum>> polynomials_m;
    std::size_t xDerivatives_m;
    std::size_t sDerivatives_m;
};

inline
    std::size_t DifferentialOperatorTwo::getXDerivatives() const {
        return xDerivatives_m;
}
inline
    std::size_t DifferentialOperatorTwo::getSDerivatives() const {
        return sDerivatives_m;
}

}

#endif
