#ifndef CLASSIC_TransportMap_HH
#define CLASSIC_TransportMap_HH

// ------------------------------------------------------------------------
// $RCSfile: TransportMap.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: TransportMap<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <iosfwd>

template <class T, int M, int N> class FMatrix;
template <class T, int N> class FVector;
template <class T, int N> class FVps;
template <class T, int N> class TransportFun;


// Template class TransportMap<T,N>
// ------------------------------------------------------------------------
/// Transport map with values of type [b]T[/b] in [b]N[/b] variables.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons
//  they are not implemented.

template <class T, int N>
class TransportMap {

public:

    /// Default constructor.
    //  Construct identity map.
    TransportMap();

    /// Convert.
    explicit TransportMap(const FVps<T, N> &rhs);

    /// Convert from matrix.
    //  The constant part is set to zero.
    //  The Transport part is filled from [b]M[/b].
    explicit TransportMap(const FMatrix<T, N, N> &M);

    /// Convert from vector.
    //  The constant part is filled from [b]V[/b].
    //  The Transport part is set to the identity.
    explicit TransportMap(const FVector<T, N> &V);

    /// Get component.
    //  Return value of component [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    const TransportFun<T, N> &getComponent(int n) const;

    /// Set component.
    //  Assign value of component [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    void setComponent(int, const TransportFun<T, N> &);

    /// Get component.
    //  Return reference to component [b]n[/b].
    //  Result is undefined, if index is out of range.
    TransportFun<T, N> &operator[](int);

    /// Get Component.
    //  Return constant reference to component [b]n[/b].
    //  Result is undefined, if index is out of range.
    const TransportFun<T, N> &operator[](int) const;

    /// Unary plus.
    TransportMap operator+() const;

    /// Unary minus.
    TransportMap operator-() const;

    /// Multiply and assign.
    TransportMap &operator*=(const TransportFun<T, N> &rhs);

    /// Divide and assign.
    //  Throw DivideError if constant part of [b]rhs[/b] is zero.
    TransportMap &operator/=(const TransportFun<T, N> &rhs);

    /// Multiply and assign.
    TransportMap &operator*=(const T &rhs);

    /// Divide and assign.
    //  Throw DivideError if [b]rhs[/b] is zero.
    TransportMap &operator/=(const T &rhs);

    /// Add.
    TransportMap &operator+=(const TransportMap &rhs);

    /// Subtract.
    TransportMap &operator-=(const TransportMap &rhs);

    /// Add and assign.
    TransportMap &operator+=(const FVector<T, N> &);

    /// Subtract and assign.
    TransportMap &operator-=(const FVector<T, N> &);

    /// Inverse.
    TransportMap inverse() const;

    /// Set to identity.
    void identity();

    /// Evaluate map at point [b]P[/b].
    FVector<T, N> constantTerm(const FVector<T, N> &P) const;

    /// Evaluate map at origin.
    /// This is equivalent to extracting constant part.
    FVector<T, N> constantTerm() const;

    /// Extract Transport terms at origin.
    //  This is equivalent to extracting Transport part.
    FMatrix<T, N, N> linearTerms() const;

    /// Substitute matrix into map.
    TransportMap substitute(const FMatrix<T, N, N> &rhs) const;

    /// Substitute map into map.
    TransportMap substitute(const TransportMap &rhs) const;

    /// Substitute map into matrix.
    TransportMap substituteInto(const FMatrix<T, N, N> &lhs) const;

protected:

    /// Size of a component.
    static const int SIZE = (N + 1) *(N + 2) / 2;

    // Representation of the TransportMap.
    TransportFun<T, N> data[N];
};


// Global functions.
// ------------------------------------------------------------------------

/// Multiply.
template <class T, int N>
TransportMap<T, N>
operator*(const TransportMap<T, N> &lhs, const TransportFun<T, N> &rhs);

/// Multiply.
template <class T, int N>
TransportMap<T, N>
operator*(const TransportFun<T, N> &lhs, const TransportMap<T, N> &rhs);

/// Multiply.
template <class T, int N>
TransportMap<T, N>
operator*(const TransportMap<T, N> &lhs, const T &rhs);

/// Multiply.
template <class T, int N>
TransportMap<T, N>
operator*(const T &lhs, const TransportMap<T, N> &rhs);

/// Multiply.
template <class T, int N>
TransportMap<T, N>
operator*(const FMatrix<T, N, N> &lhs, const TransportMap<T, N> &rhs);

/// Divide.
//  Throw DivideError, if constant part of [b]rhs[/b] is zero.
template <class T, int N>
TransportMap<T, N>
operator/(const TransportMap<T, N> &lhs, const TransportFun<T, N> &rhs);

/// Divide.
//  Throw DivideError, if [b]rhs[/b] is zero.
template <class T, int N>
TransportMap<T, N>
operator/(const TransportMap<T, N> &lhs, const T &rhs);

/// Add.
template <class T, int N>
TransportMap<T, N>
operator+(const TransportMap<T, N> &lhs, const TransportMap<T, N> &rhs);

/// Subtract.
template <class T, int N>
TransportMap<T, N>
operator-(const TransportMap<T, N> &lhs, const TransportMap<T, N> &rhs);

/// Add.
template <class T, int N>
TransportMap<T, N>
operator+(const TransportMap<T, N> &lhs, const FVector<T, N> &rhs);

/// Subtract.
template <class T, int N>
TransportMap<T, N>
operator-(const TransportMap<T, N> &lhs, const FVector<T, N> &rhs);

/// Add.
template <class T, int N>
TransportMap<T, N>
operator+(const FVector<T, N> &lhs, const TransportMap<T, N> &rhs);

/// Subtract.
template <class T, int N>
TransportMap<T, N>
operator-(const FVector<T, N> &lhs, const TransportMap<T, N> &rhs);


// Implementation.
#include "FixedAlgebra/TransportMap.hpp"

#endif // CLASSIC_TransportMap_HH
