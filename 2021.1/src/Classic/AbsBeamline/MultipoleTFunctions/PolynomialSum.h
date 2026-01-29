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

#ifndef POLYNOMIAL_SUM_H
#define POLYNOMIAL_SUM_H

/** ---------------------------------------------------------------------
  *
  * PolynomialSum describes a sum of TwoPolynomial objects.\n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  *  The polynomial of two variables @f[p(x) = a_{00} + a_{10}x + a_{11}xS(s)
  *  + ... + a_{nm}x^nS(s)^m @f] cannot be summed with other polynomials \n
  *  unless all the powers of S(s)-derivatives are identical. \n
  *  Instead all terms are just stored seperately in a list.
  */

#include <vector>
#include "TwoPolynomial.h"

namespace polynomial {

class PolynomialSum {
public:
    /** Default constructor, initialises empty sum */
    PolynomialSum();
    /** Constructor, initialises sum with a single polynomial from input
     *  \param polynomial -> The polynomial sum is initialised with this \n
     *  polynomial as the first and only term
     */
    explicit PolynomialSum(const TwoPolynomial &polynomial);
    /** Copy constructor */
    PolynomialSum(const PolynomialSum &polynomialSum);
    /** Desctructor, does nothing */
    ~PolynomialSum();
    /** Assigment operator */
    PolynomialSum& operator= (const PolynomialSum &sum);
    /** Differentiate each term wrt x */
    void differentiateX();
    /** Differentiate each term wrt s */
    void differentiateS();
    /** Multiply term with input polynomial
     *  \param poly -> Polynomial to be multiplied with this polynomial
     */
    void multiplyPolynomial(const TwoPolynomial &poly);
    /** Print polynomial, for internal debugging */
    void printPolynomial() const;
    /** Returns number of terms in the sum */
    std::size_t numberOfTerms() const;
    /** Check if term p is a zero polynomial \n
     *  Returns true if p is negative or outside the range of list range
     *  \param p -> Term index, starting with term 0
     */
    bool isPolynomialZero(const std::size_t &p) const;
    /** Truncate series in x at truncateOrder
     *  \param truncateOrder -> Highest power of x after truncation
     */
    void truncate(const std::size_t &truncateOrder);
    /** Evaluate polynomial in term p at the point (x, s) \n
     *  If p is outside list range zero is returned
     *  \param p -> Term index, starting with term 0
     *  \param x -> Point x where polynomial is evaluated
     *  \param s -> Point s where polynomial is evaluated
     */
    double evaluatePolynomial(const std::size_t &p,
                              const double &x,
                              const double &s) const;
    /** Put together all terms in the PolynomialSum by evaluating the \n
     *  S(s)-derivatives and multiply them with the polynomial coefficients \n
     *  Polynomial coefficients are now of type double 
     *  \param dSvalues -> List of evaluated S(s)-derivates, starting
     *  with zeroth (no) derivative
     *  \param finalPolynomial -> Resulting polynomial, must pass an
     *  empty vector of vectors into this function
     */
    void putSumTogether(const std::vector<double> &dSvalues,
                        std::vector<std::vector<double>> &finalPolynomial) const;
    /** Evaluate polynomial after putting the sum together into one polynomial
     *  \param x -> Point x where polynomial is evaluated
     *  \param s -> Point s where polynomial is evaluated
     *  \param dSvalues -> List of evaluated S(s)-derivatives, starting 
     *  with zeroth (no) derivative
     */
    double evaluatePolynomial2(const double &x,
                               const double &s,
                               const std::vector<double> &dSvalues) const;
    /** Add poly to the sum by concatenating the lists
     *  \param poly -> Polynomial that is added to the list
     */
    void addPolynomial(const PolynomialSum &poly);
    /** Returns lists of S(s)-derivatives in term p \n
     *  Returns empty list if p is negative or outside list range
     *  \param p -> Term index, starting with term 0
     */
    std::vector<std::size_t> getdSfactors(const std::size_t &p) const;
    /** Sort polynomialSum_m such that the TwoPolynomial objects with fewest \n
     *  S(s)-derivatives come first, in ascending order \n
     *  If any TwoPolynomial objects have identical S(s)-derivatives these 
     *  are put together by adding the polynomials
     */
    void sortTerms();
private:
    std::vector<TwoPolynomial> polynomialSum_m;
};

inline
    std::size_t PolynomialSum::numberOfTerms() const {
        return polynomialSum_m.size();
}
inline
    double PolynomialSum::evaluatePolynomial(const std::size_t &p,
                                             const double &x,
                                             const double &s) const {
        if (p >= polynomialSum_m.size()) {
            return 0.0;
        }
        return polynomialSum_m[p].evaluatePolynomial(x, s);
}

}

#endif
