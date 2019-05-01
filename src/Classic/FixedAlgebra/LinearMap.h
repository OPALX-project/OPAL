#ifndef CLASSIC_LinearMap_HH
#define CLASSIC_LinearMap_HH

// ------------------------------------------------------------------------
// $RCSfile: LinearMap.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: LinearMap<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:36 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <iosfwd>

template <class T, int M, int N> class FMatrix;
template <class T, int N> class LinearFun;
template <class T, int N> class FVector;
template <class T, int N> class FVps;


// Template class LinearMap<T,N>
// ------------------------------------------------------------------------
/// Linear map with values of type [b]T[/b] in [b]N[/b] variables.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler do the correct thing, and are not defined for speed.

template <class T, int N>
class LinearMap {

public:

    /// Default constructor.
    //  Construct identity map.
    LinearMap();

    /// Convert from general map.
    explicit LinearMap(const FVps<T, N> &rhs);

    /// Convert from matrix.
    //  The constant part is set to zero.
    //  The linear part is filled from [b]M[/b].
    explicit LinearMap(const FMatrix<T, N, N> &M);

    /// Convert from vector.
    //  The constant part is filled from [b]V[/b].
    //  The linear part is set to the identity.
    explicit LinearMap(const FVector<T, N> &V);

    /// Get component.
    //  Return value of component [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    const LinearFun<T, N> &getComponent(int n) const;

    /// Set component.
    //  Assign value of component [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    void setComponent(int, const LinearFun<T, N> &);

    /// Get component.
    //  Return reference to component [b]n[/b].
    //  Result is undefined, if index is out of range.
    LinearFun<T, N> &operator[](int);

    /// Get Component.
    //  Return constant reference to component [b]n[/b].
    //  Result is undefined, if index is out of range.
    const LinearFun<T, N> &operator[](int) const;

    /// Unary plus.
    LinearMap operator+() const;

    /// Unary minus.
    LinearMap operator-() const;

    /// Multiply and assign.
    LinearMap &operator*=(const LinearFun<T, N> &rhs);

    /// Divide and assign.
    //  Throw DivideError if constant part of [b]rhs[/b] is zero.
    LinearMap &operator/=(const LinearFun<T, N> &rhs);

    /// Multiply and assign.
    LinearMap &operator*=(const T &rhs);

    /// Divide and assign.
    //  Throw DivideError if [b]rhs[/b] is zero.
    LinearMap &operator/=(const T &rhs);

    /// Add.
    LinearMap &operator+=(const LinearMap &rhs);

    /// Subtract.
    LinearMap &operator-=(const LinearMap &rhs);

    /// Add and assign.
    LinearMap &operator+=(const FVector<T, N> &);

    /// Subtract and assign.
    LinearMap &operator-=(const FVector<T, N> &);

    /// Get a LinearMap from stream [b]is[/b].
    std::istream &get(std::istream &is);

    /// Put a LinearMap to stream [b]os[/b].
    std::ostream &put(std::ostream &os) const;

    /// Inverse.
    LinearMap inverse() const;

    /// Set to identity.
    void identity();

    /// Evaluate map at point [b]P[/b].
    FVector<T, N> constantTerm(const FVector<T, N> &P) const;

    /// Evaluate map at origin.
    /// This is equivalent to extracting constant part.
    FVector<T, N> constantTerm() const;

    /// Extract linear terms at origin.
    //  This is equivalent to extracting linear part.
    FMatrix<T, N, N> linearTerms() const;

    /// Substitute matrix into map.
    LinearMap substitute(const FMatrix<T, N, N> &rhs) const;

    /// Substitute map into map.
    LinearMap substitute(const LinearMap &rhs) const;

    /// Substitute map into matrix.
    LinearMap substituteInto(const FMatrix<T, N, N> &lhs) const;

protected:

    // Representation of the LinearMap.
    LinearFun<T, N> data[N];
};


// Global functions.
// ------------------------------------------------------------------------

/// Multiply.
template <class T, int N>
LinearMap<T, N> operator*(const LinearMap<T, N> &lhs, const LinearFun<T, N> &rhs);

/// Multiply.
template <class T, int N>
LinearMap<T, N> operator*(const LinearFun<T, N> &lhs, const LinearMap<T, N> &rhs);

/// Multiply.
template <class T, int N>
LinearMap<T, N> operator*(const LinearMap<T, N> &lhs, const T &rhs);

/// Multiply.
template <class T, int N>
LinearMap<T, N> operator*(const T &lhs, const LinearMap<T, N> &rhs);

/// Multiply.
template <class T, int N>
LinearMap<T, N> operator*(const FMatrix<T, N, N> &lhs, const LinearMap<T, N> &rhs);

/// Divide.
//  Throw DivideError, if constant part of [b]rhs[/b] is zero.
template <class T, int N>
LinearMap<T, N> operator/(const LinearMap<T, N> &lhs, const LinearFun<T, N> &rhs);

/// Divide.
//  Throw DivideError, if [b]rhs[/b] is zero.
template <class T, int N>
LinearMap<T, N> operator/(const LinearMap<T, N> &lhs, const T &rhs);

/// Add.
template <class T, int N>
LinearMap<T, N> operator+(const LinearMap<T, N> &lhs, const LinearMap<T, N> &rhs);

/// Subtract.
template <class T, int N>
LinearMap<T, N> operator-(const LinearMap<T, N> &lhs, const LinearMap<T, N> &rhs);

/// Add.
template <class T, int N>
LinearMap<T, N> operator+(const LinearMap<T, N> &lhs, const FVector<T, N> &rhs);

/// Subtract.
template <class T, int N>
LinearMap<T, N> operator-(const LinearMap<T, N> &lhs, const FVector<T, N> &rhs);

/// Add.
template <class T, int N>
LinearMap<T, N> operator+(const FVector<T, N> &lhs, const LinearMap<T, N> &rhs);

/// Subtract.
template <class T, int N>
LinearMap<T, N> operator-(const FVector<T, N> &lhs, const LinearMap<T, N> &rhs);

/// Extract LinearMap from stream [b]is[/b].
template <class T, int N>
std::istream &operator>>(std::istream &is, LinearMap<T, N> &);

/// Insert LinearMap to stream [b]os[/b].
template <class T, int N>
std::ostream &operator<<(std::ostream &os, const LinearMap<T, N> &vps);


// Implementation.
#include "FixedAlgebra/LinearMap.hpp"

#endif // CLASSIC_LinearMap_HH
