#ifndef CLASSIC_FLieGenerator_HH
#define CLASSIC_FLieGenerator_HH

// ------------------------------------------------------------------------
// $RCSfile: FLieGenerator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FLieGenerator<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2002/03/25 20:44:16 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include <complex>

// template <class T> class complex;
template <class T, int N> class FArray1D;
template <class T, int N, int M> class FMatrix;
template <class T, int N> class FTps;


// Template class FLieGenerator<T,N>.
// ------------------------------------------------------------------------
/// A representation for a homogeneous polynomial, used as a Lie generator.
//  The number of variables must be even, it is equal to 2*N, where N is
//  a template parameter.  The indices have the same range as the
//  coefficients of the same order in an FTps<double,2*N> object.

template <class T, int N>
class FLieGenerator {

public:

    /// Construct a zero generator of given order.
    explicit FLieGenerator(int);

    /// Construct a zero generator of given order by extraction
    //  from an FTps<T,2*N> object.
    FLieGenerator(const FTps<T, 2 * N> &, int);

    FLieGenerator();
    FLieGenerator(const FLieGenerator &);
    ~FLieGenerator();
    const FLieGenerator &operator=(const FLieGenerator &);
    bool operator==(const FLieGenerator &) const;


    /// Get pointer to beginning of generator.
    //  Version for non-constant generator.
    inline T *begin();

    /// Get pointer past end of generator.
    //  Version for non-constant generator.
    inline T *end();

    /// Get pointer to beginning of generator.
    //  Version for constant generator.
    inline const T *begin() const;

    /// Get pointer past end of generator.
    //  Version for constant generator.
    inline const T *end() const;

    /// Get element.
    //  Return a reference to element [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    T &operator[](int n);

    /// Get element.
    //  Return a reference to element [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    const T &operator[](int n) const;


    /// Change sign of generator.
    FLieGenerator operator-() const;

    /// Multiply by scalar and assign.
    FLieGenerator &operator*=(const T &);

    /// Divide by scalar and assign.
    FLieGenerator &operator/=(const T &);

    /// Add vector and assign.
    FLieGenerator &operator+=(const FLieGenerator &);

    /// Subtract vector and assign.
    FLieGenerator &operator-=(const FLieGenerator &);


    /// Clear all coefficients.
    void clear();

    /// Partial derivative.
    //  Return partial derivative with respect to variable [b]var[/b].
    FLieGenerator derivative(int var) const;

    /// Test for zero.
    bool isZero() const;

    /// Scale monomial-wise.
    FLieGenerator scale(const FLieGenerator &) const;

    /// Substitute matrix in Lie generator.
    //  The coefficients of the source generator have type T, the
    //  elements of the matrix have type U, and the result of multiplying
    //  a T by a U must return a U.
    template <class U>
    FLieGenerator<U, N> transform(const FMatrix<U, 2 * N, 2 * N> &) const;


    /// Return order of this generator.
    inline int getOrder() const;

    /// Return bottom index of this generator.
    inline int getBottomIndex() const;

    /// Return top index of this generator.
    inline int getTopIndex() const;

private:

    // Return size of the generator.
    static inline int getSize(int order);

    // Return bottom index for given order.
    static inline int getBottomIndex(int order);

    // Return top index for given order.
    static inline int getTopIndex(int order);

    // The order of the generator, defines its size.
    int itsOrder;

    // The bottom Index.
    int bottomIndex;

    // The top index.
    int topIndex;

    // The representation of the coefficients.
    Array1D<T> itsCoeffs;
};


// Global functions.
// ------------------------------------------------------------------------

/// Add.
template <class T, int N>
FLieGenerator<T, N>
operator+(const FLieGenerator<T, N> &, const FLieGenerator<T, N> &);

/// Subtract.
template <class T, int N>
FLieGenerator<T, N>
operator-(const FLieGenerator<T, N> &, const FLieGenerator<T, N> &);

/// Multiply by scalar.
template <class T, int N>
FLieGenerator<T, N>
operator*(const FLieGenerator<T, N> &, const T &);

/// Multiply by scalar.
template <class T, int N>
FLieGenerator<T, N>
operator*(const T &, const FLieGenerator<T, N> &);

/// Multiply by Lie generator.
template <class T, int N>
FLieGenerator<T, N>
operator*(const FLieGenerator<T, N> &, const FLieGenerator<T, N> &);

/// Divide by scalar.
template <class T, int N>
FLieGenerator<T, N>
operator/(const FLieGenerator<T, N> &, const T &);

/// Take real part of a complex generator.
template <class T, int N>
FLieGenerator<T, N>
real(const FLieGenerator<std::complex<T>, N> &);

/// Take imaginary part of a complex generator.
template <class T, int N>
FLieGenerator<T, N>
imag(const FLieGenerator<std::complex<T>, N> &);

/// Convert real generator to complex.
template <class T, int N>
FLieGenerator<std::complex<T>, N>
toComplex(const FLieGenerator<T, N> &);

/// Poisson bracket of two Lie generators.
template <class T, int N> FLieGenerator<T, N>
PoissonBracket(const FLieGenerator<T, N> &, const FLieGenerator<T, N> &);

/// Output.
template <class T, int N>
std::ostream &operator<<(std::ostream &, const FLieGenerator<T, N> &);

#include "FixedAlgebra/FLieGenerator.hpp"

#endif // CLASSIC_FLieGenerator_HH
