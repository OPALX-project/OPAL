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

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

/** ---------------------------------------------------------------------
  *
  * Polynomial describes a polynomial of one variable.\n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  *  The polynomial @f[p(x) = a_0 + a_1x + ... + a_nx^n @f] is stored as a \n
  *  list (a_0, a_1, ..., a_n).\n
  *  BUG: For large n the integer type might overflow. If you need this many \n
  *  terms change all int to long int. \n
  */

#include <vector>

namespace polynomial {

class Polynomial {
public:
    /** Default constructor, initialises zero polynomial */
    Polynomial();
    /** Initialises polynomial with input coefficients, starting with the
     *  const term \n
     *  If coefficients is an empty vector, a zero polynomial with one
     *  constant zero term is created 
     *  \param coefficients -> List of polynomial coefficients, starting
     *  with the const term
     */
    explicit Polynomial(const std::vector<int> &coefficients);
    /** Copy constructor */
    Polynomial(const Polynomial &poly);
    /** Destructor, does nothing */
    ~Polynomial();
    /** Assigment operator */
    Polynomial& operator= (const Polynomial &poly);
    /** Differentiate polynomial using \f$\frac{d}{dx}x^n = nx^{n - 1}\f$ */
    void differentiatePolynomial();
    /** Multiply polynomial with input polynomial by multiplying
     *  term by term \n
     *  maxXorder_m is adjusted accordingly and coefficients_m is resized 
     *  \param poly -> Polynomial to be multiplied with this polynomial
     */
    void multiplyPolynomial(const Polynomial &poly);
    /** Add polynomial to input polynomial \n
     *  If input polynomial has higher order terms,
     *  then maxXorder_m is adjusted and coefficients_m is resized
     *  \param poly -> Polynomial to be added to this polynomial
     */
    void addPolynomial(const Polynomial &poly);
    /** Print polynomial for internal debugging */
    void printPolynomial() const;
    /** Returns coefficient in front of x^order
     *  \param order -> The power of x belonging to this coefficient
     */
    int getCoefficient(const std::size_t &order) const;
    /** Set the coefficient in front of x^order \n
     *  If order is larger than maxXorder_m maxXorder_m is adjusted
     *  and coefficients_m is resized
     *  \param coefficient -> Coefficient we wish to insert
     *  \param order -> The power of x belonging to this coefficient
     */
    void setCoefficient(const int &coefficient, const std::size_t &order);
    /** Returns the highest order of x, does nothing if order is negative */
    int getMaxXorder() const;
    /** Set the highest order of x to maxXorder, filling in zeros if necessary
     *  \param maxXorder -> Highest power of x
     */
    void setMaxXorder(const std::size_t &maxXorder);
    /** Set polynomial to zero */
    void setZero();
    /** Evaluate the polynomial
     *  \param x -> The point at which polynomial is evaluated
     */
    double evaluatePolynomial(const double &x) const;
    /** Truncate polynomial to truncateOrder \n
     *  If truncateOrder is greater than maxXorder_m, maxXorder_m is adjusted
     *  and coefficients is resized with zero coefficients at the end
     *  \param truncateOrder -> Highest order of x after truncation
     */
    void truncate(const std::size_t &truncateOrder);
private:
    std::size_t maxXorder_m;
    std::vector<int> coefficients_m;
};

inline
    int Polynomial::getCoefficient(const std::size_t &order) const {
        if (order > maxXorder_m) {
            return 0;
    }
        return coefficients_m[order];
}
inline
    int Polynomial::getMaxXorder() const{
        return maxXorder_m;
}
inline
    void Polynomial::setMaxXorder(const std::size_t &maxXorder) {
        coefficients_m.resize(maxXorder + 1, 0);
        maxXorder_m = maxXorder;
}
inline
    void Polynomial::setZero() {
        maxXorder_m = 0;
        coefficients_m.resize(1);
        coefficients_m[0] = 0;
}
inline
    void Polynomial::truncate(const std::size_t &truncateOrder) {
        if (truncateOrder < maxXorder_m) {
            coefficients_m.resize(truncateOrder + 1);
            maxXorder_m = truncateOrder;
        }
}

}

#endif
