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

#ifndef RECURSION_RELATION_H
#define RECURSION_RELATION_H

/** ---------------------------------------------------------------------
  *
  * RecursionRelation describes the differential operator used to find the \n 
  * coefficients in the expansion of the magnetic scalar potential \n
  * It contains member functions for extracting all information required to \n
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
  *  \Big(\frac{1}{\rho(1 &+& x/\rho)}\frac{\partial}{\partial x} 
  *  + \frac{\partial^2}{\partial x^2}
  *  &+& \frac{1}{(1 + x/\rho)^2}\frac{\partial^2}{\partial s^2}\Big)^n
  *  \f}
  *  and it can be initialised to any power of x and up to any n.
  */

#include "DifferentialOperator.h"
#include "Polynomial.h"

namespace polynomial {

class RecursionRelation {
public:
    /** Default constructor, initialises identity operator */
    RecursionRelation();
    /** Constructor, initialises the operator
     *  \f{eqnarray*} {
     *  \Big(\frac{1}{\rho(1 + x/\rho)}\frac{\partial}{\partial x}
     *  + \frac{\partial^2}{\partial x^2}
     *  + \frac{1}{(1 + x/\rho)^2}\frac{\partial^2}{\partial s^2}\Big)^n
     *  \f}
     *  where power = n. The denominators are expanded in x and \n
     *  highestXorder is the highest power of x \n
     *  \param power -> Number of times the differential operator is applied
     *  \param highestXorder -> Highest order of x before truncation
     */
    RecursionRelation(const std::size_t &power,
                      const std::size_t &highestXorder);
    /** Copy constructor */
    RecursionRelation(const RecursionRelation &doperator);
    /** Desctructor, does nothing */
    ~RecursionRelation();
    /** Assigment operator */
    RecursionRelation& operator= (const RecursionRelation &recursion);
    /** Print operator, used for internal debugging */
    void printOperator() const;
    /** Truncate series in x at highestXorder
     *  \param highestXorder -> Highest order of x after truncation
     */
    void truncate(std::size_t highestXorder);
    /** Increase power n by one */
    void applyOperator();
    /** Apply a differential operator in x */
    void differentiateX();
    /** Returns highest derivative of x */
    std::size_t getMaxXDerivatives() const;
    /** Returns highest derivative of s */
    std::size_t getMaxSDerivatives() const;
    /** Evaluates polynomial
     *  \param x -> Point at which polynomial is evaluated
     *  \param xDerivative -> Number of x-derivatives
     *  \param sDerivative -> Number of s-derivatives
     */
    double evaluatePolynomial(const double &x,
                              const std::size_t &xDerivative,
                              const std::size_t &sDerivative) const;
    /** Check if polynomial with x x-derivatives and s s-derivatives is zero
     *  \param xDerivative -> Number of x-derivatives
     *  \param sDerivative -> Number of s-derivatives
     */
    bool isPolynomialZero(const std::size_t &x, const std::size_t &s) const;
    /** Change number of x-derivatives to xDerivatives
     *  \param xDerivative -> Number of x-derivatives
     */  
    void resizeX(const std::size_t &xDerivatives);
    /** Change number of s-derivatives to sDerivatives
     *  \param sDerivative -> Number of s-derivatives
     */
    void resizeS(const std::size_t &sDerivatives);
private:
    DifferentialOperator operator_m;
    std::size_t power_m;
    std::size_t highestXorder_m;
};

inline
    void RecursionRelation::printOperator() const {
        operator_m.printOperator();
}
inline
    void RecursionRelation::differentiateX() {
        operator_m.differentiateX();
}
inline
    std::size_t RecursionRelation::getMaxXDerivatives() const {
        return operator_m.getXDerivatives();
}
inline
    std::size_t RecursionRelation::getMaxSDerivatives() const {
        return operator_m.getSDerivatives();
}
inline
    double RecursionRelation::evaluatePolynomial(
                              const double &x,
                              const std::size_t &xDerivative,
                              const std::size_t &sDerivative) const {
        return operator_m.evaluatePolynomial(x, xDerivative, sDerivative);
}
inline
    bool RecursionRelation::isPolynomialZero(const std::size_t &x,
                                             const std::size_t &s) const {
        return operator_m.isPolynomialZero(x, s);
}
inline
    void RecursionRelation::resizeX(const std::size_t &xDerivatives) {
        operator_m.resizeX(xDerivatives);
}
inline
    void RecursionRelation::resizeS(const std::size_t &sDerivatives) {
        operator_m.resizeS(sDerivatives);
}

}

#endif
