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

#ifndef TWO_POLYNOMIAL_H
#define TWO_POLYNOMIAL_H

/** ---------------------------------------------------------------------
  *
  * TwoPolynomial describes a polynomial of two variables. In addition, \n
  * TwoPolynomial also stores a list of derivatives of the fringe field S(s), \n
  * because the second variable is S(s). Every s-derivative will therefore \n
  * give various terms of S(s)-derivatives from the chain rule. \n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  *  The polynomial @f[p(x) = a_{00} + a_{10}x + a_{11}xS(s) + ...
  *  + a_{nm}x^nS(s)^m $f] is \n stored as a two-dimensional list \n
  *  (a_{00}, a_{10}, ..., a_{n0}; a_{01}, ...). \n
  *  In addition a list @f[b_k] of powers of S(s)-derivatives is stored,
  *  representing @f[\prod_k\frac{dS(s)}{ds}^{b_k}]
  *  A TwoPolynomial object therefore represents a product of a polynomial \n
  *  in x and S(s), multiplied by derivatives of S(s), and it is the powers \n
  *  that are stored in this list
  *  BUG: For large n the integer type might overflow. If you need this \n
  *  many terms change all int to long int. \n
  */

#include <vector>

namespace polynomial {

/** Struct for storing vector lengths used in multiplyPolynomial()
 *  \param nx -> Size of first dimension of resulting array
 *  \param ny -> Size of second dimension of resulting array
 *  \param nx1 -> Size of first dimension of first array
 *  \param ny1 -> Size of second dimension of first array
 *  \param nx2 -> Size of first dimension of second array
 *  \param ny2 -> Size of second dimension of second array
 */
struct vectorLengths {
    std::size_t nx;
    std::size_t ny;
    std::size_t nx1;
    std::size_t ny1;
    std::size_t nx2;
    std::size_t ny2;
};

class TwoPolynomial {
public:
    /** Default constructor, initialises zero polynomial */
    TwoPolynomial();
    /** Constructor, initialises polynomial with the input coefficients \n
     *  If either list dimensions are empty, a zero polynomial is initialised
     *  \param coefficients -> 2D list of polynomial coefficients, with
     *  x powers as first argument and s powers as second argument, starting
     *  with a constant term
     */
    explicit TwoPolynomial(const std::vector<std::vector<int>> &coefficients);
    /** Copy constructor */
    TwoPolynomial(const TwoPolynomial &poly);
    /** Assigment operator */
    TwoPolynomial& operator= (const TwoPolynomial &poly);
    /** Destructor, does nothing */
    ~TwoPolynomial();
    /** Differentiate polynomial wrt x */
    void differentiateX();
    /** Differentiate polynomial wrt s */
    void differentiateS();
    /** Multiply polynomial with input constant
     *  \param constant -> Constant to be multiplied with polynomial
     */
    void multiplyConstant (const int &constant);
    /** Multiply S(s)-derivatives part with input polynomial
     *  \param poly -> Polynomial to be multiplied with current polynomial
     */
    void multiplydSfactors(const TwoPolynomial &poly);
    /** Convert two 2D vector<int> array into 1D double C-array by indexing
     *  \param vec1 -> Pointer to the first 1D array
     *  \param vec2 -> Pointer to the second 1D array
     *  \param nn -> Vector sizes
     *  \param poly -> Polynomial to be multiplied with current polynomial
     */
    void convert2Dto1Darray(double *vec1, double *vec2,
                            const vectorLengths &nn,
                            const TwoPolynomial &poly) const;
    /** Perform a 2D convolution on two 1D indexed arrays using FFT \n
     *  First a 2D FFT is done by doing 1D FFT on each row and and column \n
     *  Convolution in Fourier space is a simple product \n
     *  Then the inverse 2D FFT is done on first list
     *  \param vec1 -> First array, also resulting array
     *  \param vec2 -> Second array
     *  \param nx -> Size of first dimension of resulting array
     *  \param ny -> Size of second dimension of resulting array
     */
    void convolution(double *vec1, double *vec2,
                     const std::size_t &nx, const std::size_t &ny) const;
    /** Convert 1D double indexed C-array to 2D vector<int> array by rounding
     *  \param vec1 -> Pointer to resulting 1D array
     *  \param nn -> Vector sizes
     */
    void convert1Dto2Darray(double *vec1, const vectorLengths &nn);
    /** Multiply polynomial with input polynomial using FFT and convolution
     *  \param poly -> Polynomial to be multiplied with current polynomial
     */
    void multiplyPolynomial(const TwoPolynomial &poly);
    /** Print polynomial, for internal debugging */
    void printPolynomial() const;
    /** Returns coefficient with power Xorder in x and Sorder in s
     *  \param Xorder -> Power of x
     *  \param Sorder -> Power of s
     */
    int getCoefficient(const std::size_t &Xorder,
                       const std::size_t &Sorder) const;
    /** Set coefficient of term with Xorder powers of x and Sorder powers of s
     *  \param coefficient -> Value assigned to polynomial coefficient
     *  \param Xorder -> Power of x
     *  \param Sorder -> Power of s
     */
    void setCoefficient(const int &coefficient,
                        const std::size_t &Xorder,
                        const std::size_t &Sorder);
    /** Returns highest power of x */
    std::size_t getMaxXorder() const;
    /** Returns highest power of s */
    std::size_t getMaxSorder() const;
    /** Set highest power of x
     *  \param maxXorder -> Highest power of x
     */
    void setMaxXorder(const std::size_t &maxXorder);
    /** Set highest power of s
     *  \param maxSorder -> Highest power of s
     */
    void setMaxSorder(const std::size_t &maxSorder);
    /** Return highest derivative of S(s) */
    std::size_t getNdSfactors() const;
    /** Returns the list of powers of S(s)-derivatives */
    std::vector<std::size_t> getdSfactors() const;
    /** Set the list of S(s)-derivatives
     *  \param dSfactors -> List of S(s)-derivatives
     */
    void setdSfactors(const std::vector<std::size_t> &dSfactors);
    /** Set polynomial to zero */
    void setZero();
    /** Check if polynomial is zero */
    bool isZero() const;
    /** Evaluate polynomial at point (x, s) 
     *  \param x -> x-coordinate where polynomial is evaluated
     *  \param s -> s-coordinate where polynomial is evaluated
     */
    double evaluatePolynomial(const double &x, const double &s) const;
    /** Truncate polynomial in x at truncateOrder 
     *  \param truncateOrder -> Highest power of x after truncation
     */
    void truncate(const std::size_t &truncateOrder);
    /** Add two polynomials \n
     *  The S(s)-derivatives are untouched, the user must ensure 
     *  that these are identical!
     *  \param poly -> TwoPolynomial added to this polynomial
     */
    void addPolynomial(const TwoPolynomial &poly);
    /** Operator < definition
     *  First the TwoPolynomial with the longest dSfactors_m list is greatest \n
     *  If dSfactors_m have the same length, each list element is compared
     *  individually, starting with the first element (lowest derivative)
     *  \param poly -> The TwoPolynomial to be compared with
     */
    friend bool operator < (const TwoPolynomial &left,
                            const TwoPolynomial &right);
    /** Operator == definition \n
     *  Returns true if dSfactors_m are identical \n
     *  Only uses the < operator
     *  \param poly -> The TwoPolynomial to be compared with
     */
    friend bool operator == (const TwoPolynomial &left,
                             const TwoPolynomial &right);
    /** PolynomialSum is a sum of TwoPolynomials, useful to have as friend */
    friend class PolynomialSum;
private:
    std::size_t maxXorder_m;
    std::size_t maxSorder_m;
    std::vector<std::vector<int>> coefficients_m;
    std::vector<std::size_t> dSfactors_m;
};

inline
    std::size_t TwoPolynomial::getMaxXorder() const {
        return maxXorder_m;
}
inline
    std::size_t TwoPolynomial::getMaxSorder() const {
        return maxSorder_m;
}
inline
    std::size_t TwoPolynomial::getNdSfactors() const {
        return dSfactors_m.size();
}
inline
    std::vector<std::size_t> TwoPolynomial::getdSfactors() const {
        return dSfactors_m;
}
inline
    void TwoPolynomial::setdSfactors(
                        const std::vector<std::size_t> &dSfactors) {
        dSfactors_m = dSfactors;
}
inline
    bool TwoPolynomial::isZero() const {
        return (maxXorder_m == 0 && maxSorder_m == 0 &&
                coefficients_m[0][0] == 0);
}
inline
    void TwoPolynomial::truncate(const std::size_t &truncateOrder) {
        if (truncateOrder < maxXorder_m) {
            coefficients_m.resize(truncateOrder + 1);
            maxXorder_m = truncateOrder;
        }
}

}

#endif
