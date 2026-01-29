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

#ifndef RECURSION_RELATION_TWO_H
#define RECURSION_RELATION_TWO_H

/** ---------------------------------------------------------------------
  *
  * RecursionRelationTwo describes the differential operator used to find \n
  * the coefficients in the expansion of the magnetic scalar potential. It \n
  * contains member functions for extracting all information required to \n
  * reconstruct the differential operator and evaluate its terms.\n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  * The operator of interest is
  *  \f{eqnarray*} {
  *  \Big(\frac{1}{\rho(s)(1 &+& x/\rho(s))}\frac{\partial}{\partial x} +
  *  \frac{\partial^2}{\partial x^2} &+& \frac{1}{1 +
  *  x/\rho(s)}\frac{\partial}{\partial s}\Big(\frac{1}{1 +
  *  x/\rho(s)}\frac{\partial}{\partial s}\Big)\Big)^n
  *  \f}
  *  and it can be initialised to any power of x and up to any n.
  */

#include <vector>
#include "DifferentialOperatorTwo.h"
#include "PolynomialSum.h"
#include "TwoPolynomial.h"

namespace polynomial {

class RecursionRelationTwo {
public:
    /** Default constructor, initialises identity operator */
    RecursionRelationTwo();
    /** Constructor, initialises the operator
     *  \f{eqnarray*} {
     *  \Big(\frac{1}{\rho(s)(1 + x/\rho(s))}\frac{\partial}{\partial x} +
     *  \frac{\partial^2}{\partial x^2} + \frac{1}{1 + x/\rho(s)}
     *  \frac{\partial}{\partial s}\Big(\frac{1}{1 + x/\rho(s)}
     *  \frac{\partial}{\partial s}\Big)\Big)^n
     *  \f}
     *  where power = n. The denominators are expanded in x and highestXorder \n
     *  is the highest power of x \n
     *  \param power -> Value of n of the coefficient fn
     *  \param highestXorder -> Highest power of x
     */
    RecursionRelationTwo(const std::size_t &power,
                         const std::size_t &highestXorder);
    /** Copy constructor */
    RecursionRelationTwo(const RecursionRelationTwo &doperator);
    /** Desctructor, does nothing */
    ~RecursionRelationTwo();
    /** Assigment operator */
    RecursionRelationTwo& operator= (const RecursionRelationTwo &recursion);
    /** Print operator, used for internal debugging */
    void printOperator() const;
    /** Truncate series in x at highestXorder
     *  \param highestXorder -> Highest power of x after truncation
     */
    void truncate(std::size_t highestXorder);
    /** Applies another differential operator to the existing operator */
    void applyOperator();
    /** Apply a differential operator in x */
    void differentiateX();
    /** Apply a differential operator in s */
    void differentiateS();
    /** Returns highest derivative of x */
    std::size_t getMaxXDerivatives() const;
    /** Returns highest derivative of s */
    std::size_t getMaxSDerivatives() const;
    /** Evaluates polynomial term p with xDerivative x-derivatives \n
     *  and sDerivative s-derivatives
     *  If xDerivative, sDerivative or term are out or range
     *  then zero is returned \n
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
    /** Check if polynomial term p with x x-derivatives
     *  and s s-derivatives is zero \n
     *  If x, s or term are negative true is returned
     *  \param x -> Number of x-derivatives
     *  \param s -> Number of s-derivatives
     *  \param term -> Term index
     */
    bool isPolynomialZero(const std::size_t &x,
                          const std::size_t &s,
                          const std::size_t &term) const;
    /** Change number of x-derivatives to xDerivatives
     *  \param xDerivatives -> Number of x-derivatives after resize
     */
    void resizeX(const std::size_t &xDerivatives);
    /** Change number of s-derivatives to sDerivatives
     *  \param sDerivatives -> Number of s-derivatives after resize
     */
    void resizeS(const std::size_t &sDerivatives);
    /** Returns number of terms in the sum with xDerivative x-derivatives
     *  and sDerivative s-derivatives \n
     *  If xDerivative or sDerivative is negative zero is returned
     *  \param x -> Number of x-derivatives
     *  \param s -> Number of s-derivatives
     */
    std::size_t numberOfTerms(const std::size_t &xDerivative,
                              const std::size_t &sDerivative) const;
    /** Returns list of S(s)-derivatives in term p with xDerivative
     *  x-derivatives and sDerivative s-derivatives \n
     *  If xDerivative or sDerivative are out of range
     *  an empty list is returned \n
     *  \param x -> Number of x-derivatives
     *  \param s -> Number of s-derivatives
     *  \param term -> Term index
     */
    std::vector<std::size_t> getdSfactors(const std::size_t &xDerivative,
                                  const std::size_t &sDerivative,
                                  const std::size_t &p) const;
    /** Sort the terms in each sum, with fewest S(s)-derivatives first,
     *  then in increasing powers
     */
    void sortTerms();
private:
    DifferentialOperatorTwo operator_m;
    std::size_t power_m;
    std::size_t highestXorder_m;
};

/*inline
    void RecursionRelationTwo::printOperator() const {
        operator_m.printOperator();
}*/
inline
    void RecursionRelationTwo::truncate(std::size_t highestXorder) {
        highestXorder_m = highestXorder;
        operator_m.truncate(highestXorder);
}
inline
    void RecursionRelationTwo::differentiateX() {
        operator_m.differentiateX();
}
inline
    void RecursionRelationTwo::differentiateS() {
        operator_m.differentiateS();
}
inline
    std::size_t RecursionRelationTwo::getMaxXDerivatives() const {
        return operator_m.getXDerivatives();
}
inline
    std::size_t RecursionRelationTwo::getMaxSDerivatives() const {
        return operator_m.getSDerivatives();
}
inline
    double RecursionRelationTwo::evaluatePolynomial(
                                 const double &x,
                                 const double &s,
                                 const std::size_t &xDerivative,
                                 const std::size_t &sDerivative,
                                 const std::vector<double> &dSvalues) const {
        return operator_m.evaluatePolynomial(x, s,
                                             xDerivative, sDerivative,
                                             dSvalues);
}
inline
    bool RecursionRelationTwo::isPolynomialZero(const std::size_t &x,
                                                const std::size_t &s,
                                                const std::size_t &term) const {
        return operator_m.isPolynomialZero(x, s, term);
}
inline
    void RecursionRelationTwo::resizeX(const std::size_t &xDerivatives) {
        operator_m.resizeX(xDerivatives);
}
inline
    void RecursionRelationTwo::resizeS(const std::size_t &sDerivatives) {
        operator_m.resizeS(sDerivatives);
}
inline
    std::size_t RecursionRelationTwo::numberOfTerms(
                                      const std::size_t &xDerivative,
                                      const std::size_t &sDerivative) const {
        return operator_m.numberOfTerms(xDerivative, sDerivative);
}
inline
    std::vector<std::size_t> RecursionRelationTwo::getdSfactors(
                                           const std::size_t &xDerivative,
                                           const std::size_t &sDerivative,
                                           const std::size_t &p) const {
        return operator_m.getdSFactors(xDerivative, sDerivative, p);
}
inline
    void RecursionRelationTwo::sortTerms() {
        operator_m.sortTerms();
}

}

#endif
