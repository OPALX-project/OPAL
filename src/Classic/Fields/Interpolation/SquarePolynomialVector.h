/* 
 *  Copyright (c) 2015, Chris Rogers
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

#ifndef SquarePolynomialVector_hh
#define SquarePolynomialVector_hh 1

#include "Fields/Interpolation/MMatrix.h"
#include "Fields/Interpolation/MVector.h"
#include "Fields/Interpolation/PolynomialCoefficient.h"

#include <map>

namespace interpolation {

/** SquarePolynomialVector, an arbitrary order polynomial vector class.
 *  \brief SquarePolynomialVector describes a vector of multivariate polynomials
 *
 *  Consider the vector of multivariate polynomials
 *  \f$y_i = a_0 + Sum (a_j x^j)\f$
 *  i.e. maps a vector \f$\vec{x}\f$ onto a vector \f$\vec{y}\f$ with\n
 *  \f$\vec{y} = a_0 + sum( a_{j_1j_2...j_n} x_1^{j_1} x_2^{j_2} ... x_n^{j_n})\f$.
 *  Also algorithms to map a single index J into \f$j_1...j_n\f$.\n
 *  \n
 *  PolynomialVector represents the polynomial coefficients as a matrix of
 *  doubles so that\n
 *  \f$ \vec{y} = \mathbf{A} \vec{w} \f$\n
 *  where \f$\vec{w}\f$ is a vector with elements given by \n
 *  \f$ w_i = x_1^{i_1} x_2^{i_2} ... x_n^{i_n} \f$ \n
 *  So the index \f$ i \f$ is actually itself a vector. The vectorisation of 
 *  \f$ i \f$ is handled by IndexByPower and IndexByVector methods.
 *  IndexByPower gives an index like:\n
 *  \f$ w_i = x_1^{i_1} x_2^{i_2} ... x_n^{i_n} \f$ \n
 *  While IndexByVector gives an index like:\n
 *  \f$ w_i = x_{i_1}x_{i_2} \ldots x_{i_n} \f$ \n
 *  \n
 *
 *  Nb: it is a SquarePolynomialVector because coefficients include all
 *  polynomial coefficients with \f$ i_ j <= n \f$; coefficients sit within an
 *  n-dimensional square. The distinction should be made with a PolynomialVector
 *  where coefficients include all polynomial coefficients with 
 *  \f$ \Sigma(i_j) <= n. \f$
 * 
 *  Nb2: It is quite tricky to get right - mostly an excercise in indexing, but
 *  that is very easy to get wrong.
 */

class SquarePolynomialVector {
  public:
    /** Default constructor
     *
     *  Leaves everything empty, call SetCoefficients to set up properly.
     */
    SquarePolynomialVector();

    /** Copy constructor
     */
    SquarePolynomialVector (const SquarePolynomialVector& pv);

    SquarePolynomialVector& operator=(const SquarePolynomialVector&) = default;

    /** Construct polynomial vector passing polynomial coefficients as a matrix.
     */
    SquarePolynomialVector (int pointDimension,
                            MMatrix<double> polynomialCoefficients);

    /** Construct polynomial vector passing polynomial coefficients as a list of
     *  PolynomialCoefficient objects.
     *
     *  Any coefficients missing from the vector are set to 0. Maximum order of
     *  the polynomial is given by the maximum order of the coefficients in the
     *  vector.
     */
    SquarePolynomialVector (std::vector<PolynomialCoefficient> coefficients);

    /** Destructor - no memory allocated so doesn't do anything */
    ~SquarePolynomialVector() {;}

    /** Reinitialise the polynomial vector with new point (x) dimension and
     *  coefficients
     */
    void             SetCoefficients(int pointDim, MMatrix<double> coeff);

    /** Reinitialise the polynomial vector with new point (x) dimension and
     *  coefficients.
     *
     *  Any coefficients missing from the vector are set to 0. Maximum order of
     *  the polynomial is given by the maximum order of the coefficients in the
     *  vector.
     */
    void             SetCoefficients(std::vector<PolynomialCoefficient> coeff);

    /** Set the coefficients.
     *
     *  This method can't be used to change point dimension or value dimension -
     *  any range mismatch will be ignored.
     */
    void             SetCoefficients(MMatrix<double> coeff);

    /** Return the coefficients as a matrix of doubles. */
    MMatrix<double>                    GetCoefficientsAsMatrix() const
    {return _polyCoeffs;}

    /** Fill value with \f$ y_i \f$ at some set of \f$ x_i \f$ (point).
     *
     *  point should be array of length PointDimension().
     *  value should be array of length ValueDimension().
     */
    void  F(const double*   point,    double* value) const;

    /** Fill value with \f$ y_i \f$ at some set of \f$ x_i \f$ (point).
     *
     *  point should be vector of length PointDimension().
     *  value should be vector of length ValueDimension().
     *  Note that there is no bound checking here.
     */
    void  F(const MVector<double>& point, MVector<double>& value) const;

    /** Length of the input point (x) vector. */
    unsigned int PointDimension()         const {return _pointDim;}

    /** Length of the output value (y) vector. */
    unsigned int ValueDimension() const {return _polyCoeffs.num_row();}

    /** Index of highest power - 0 is const, 1 is linear, 2 is quadratic etc...
     */
    unsigned int PolynomialOrder() const;

    /** Polymorphic copy constructor.
     *
     *  This is a special copy constructor for inheritance structures
     */
    SquarePolynomialVector* Clone() const
    {return new SquarePolynomialVector(*this);}
 
    /** Make a vector like \f$(c, x, x^2, x^3...)\f$.
     *  polyVector should be of size NumberOfPolynomialCoefficients().
     *  could be static but faster as member function (use lookup table for
     *  _polyKey).
     */
    MVector<double>& MakePolyVector(const MVector<double>& point,
                                    MVector<double>& polyVector) const;

    /** Make a vector like \f$(c, x, x^2, x^3...)\f$.
     *
     *  PolyVector should be of size NumberOfPolynomialCoefficients().
     *  Could be static but faster as member function (use lookup table for
     *  _polyKey).
     */
    double* MakePolyVector(const double* point, double* polyVector) const;

    /** Transforms from a 1d index of polynomial coefficients to an nd index.
     *
     *  This is slow - you should use it to build a lookup table.
     *  For polynomial term \f$x_1^i x_2^j ... x_d^n\f$ index like [i][j] ... [n]
     *  e.g. \f$x_1^4 x_2^3 = x_1^4 x_2^3 x_3^0 x_4^0\f$ = {4,3,0,0}.
     */
    static std::vector<int> IndexByPower (int index, int nInputVariables);

    /** Transforms from a 1d index of polynomial coefficients to an nd index.
     *
     *  This is slow - you should use it to build a lookup table.
     *  For polynomial term \f$x_i x_j...x_n\f$ index like [i][j] ... [n]
     *  e.g. \f$x_1^4 x_2^3\f$ = {1,1,1,1,2,2,2}
     */
    static std::vector<int> IndexByVector(int index, int nInputVariables);

    /** Returns the number of coefficients required for an arbitrary dimension,
     *  order polynomial
     *  e.g. \f$a_0 + a_1 x + a_2 y + a_3 z + a_4 x^2 + a_5 xy + a_6 y^2 + a_7 xz + a_8 yz + a_9 z^2\f$
     *  => NumberOfPolynomialCoefficients(3,2) = 9 (indexing starts from 0).
     */
    static unsigned int NumberOfPolynomialCoefficients(int pointDimension, int order);

    /** Write out the PolynomialVector (effectively just polyCoeffs). */
    friend std::ostream& operator<<(std::ostream&,  const SquarePolynomialVector&);

    /** Control the formatting of the polynomial vector.
     *
     *  If PrintHeaders is true, then every time I write a PolynomialVector it
     *  will write the header also (default).
     */
    static void PrintHeaders(bool willPrintHeaders) {_printHeaders = willPrintHeaders;}

    /** Write out the header (PolynomialByPower index for each element). */
    void PrintHeader(std::ostream& out,
                     char int_separator='.',
                     char str_separator=' ',
                     int length=14,
                     bool pad_at_start=true) const;

    /** Print a sequence: with some string that separates elements of index and
     *  some character that pads
     *
     *  The total vector length is the total length of the output, set to < 1 to
     *  ignore (means no padding)
     *  pad_at_start determines whether to pad at start (i.e. right align) or 
     *  end (i.e. left align)
     */
    template < class Container >
    static void PrintContainer(std::ostream& out,
                               const Container& container,
                               char T_separator,
                               char str_separator,
                               int length,
                               bool pad_at_start);

private:
    static void IndexByPowerRecursive(
                                std::vector<int> check,
                                size_t check_index,
                                size_t poly_power,
                                std::vector<std::vector<int> >& nearby_points);

    int                                _pointDim;
    MMatrix<double>                    _polyCoeffs;
    static std::vector< std::vector< std::vector<int> > > _polyKeyByPower;
    static std::vector< std::vector< std::vector<int> > > _polyKeyByVector;
    static bool                        _printHeaders;
};

std::ostream& operator<<(std::ostream&, const SquarePolynomialVector&);

}
#endif // SquarePolynomialVector_hh




