#ifndef CLASSIC_FVps_HH
#define CLASSIC_FVps_HH

// ------------------------------------------------------------------------
// $RCSfile: FVps.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.6 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: FVps<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:06:20 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include <algorithm>
#include <iosfwd>

template <class T> class Array1D;
template <class T, int M, int N> class FMatrix;
template <class T, int N> class FTps;
template <class T, int N> class FVector;
template <class T, int N> class LinearMap;
template <class T, int N> class TransportMap;
template <class T, int N> class FArray1D;

// Template class FVps<T,N>
// ------------------------------------------------------------------------
/// Vector truncated power series in n variables.

template <class T, int N>
class FVps {

public:

    /// Constructor.
    //  Construct zero map with given order [b]order[/b]
    //  and truncation [b]trunc[/b].
    explicit FVps(int minOrder, int maxOrder, int trcOrder);

    /// Convert from matrix.
    //  The constant part is set to zero.
    //  The linear part is filled from [b]M[/b].
    explicit FVps(const FMatrix<T, N, N> &M);

    /// Convert from vector.
    //  The constant part is filled from [b]V[/b].
    explicit FVps(const FVector<T, N> &V);

    /// Convert from linear map.
    explicit FVps(const LinearMap<T, N> &rhs);

    /// Convert from second-order map.
    explicit FVps(const TransportMap<T, N> &rhs);

    FVps();
    FVps(const FVps &);
    ~FVps();
    const FVps &operator=(const FVps &);

    /// Set to identity.
    void identity();

    /// Set to zero.
    void zero();

    /// Get component.
    //  Return value of component [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    const FTps<T, N> &getComponent(int n) const;

    /// Set component.
    //  Assign value of component [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    void setComponent(int, const FTps<T, N> &);

    /// Get Component.
    //  Return constant reference to component [b]n[/b].
    //  Result undefined for index out of range.
    inline const FTps<T, N> &operator[](int) const;

    /// Get (Set) component.
    //  Return reference to component [b]n[/b].
    //  Result undefined for index out of range.
    inline FTps<T, N> &operator[](int);

    /// Get dimension.
    //  Return the number of components.
    inline int getDimension() const;

    /// Get number of variables.
    //  This is the same in all components.
    inline int getVariables() const;

    /// Get lowest order contained in any component.
    int getMinOrder() const;

    /// Set minimum order.
    //  If necessary, this function will insert zeroes or modify the
    //  maximum order.  It will not modify the truncation order.
    void setMinOrder(int order);

    /// Get highest order contained in any component.
    int getMaxOrder() const;

    /// Set maximum order.
    //  If necessary, this function will insert zeroes or modify the minimum
    //  or maximum orders.  It will not modify the truncation order.
    void setMaxOrder(int order);

    /// Get highest order contained in any component.
    int getTopOrder() const;

    /// Get lowest truncation order in any component.
    int getTruncOrder() const;

    /// Set truncation order for all components.
    void setTruncOrder(int order);

    /// Extract given range of orders, with truncation.
    FVps filter(int minOrder, int maxOrder, int trcOrder = (FTps<T, N>::EXACT)) const;

    /// Truncate.
    //  Limit truncation order to [b]trunc[/b].
    FVps truncate(int trunc);

    /// Unary plus.
    FVps operator+() const;

    /// Unary minus.
    FVps operator-() const;

    /// Add and assign.
    FVps &operator+=(const FVps &rhs);

    /// Subtract and assign.
    FVps &operator-=(const FVps &rhs);

    /// Add and assign.
    FVps &operator+=(const FVector<T, N> &);

    /// Subtract and assign.
    FVps &operator-=(const FVector<T, N> &);

    /// Multiply and assign.
    FVps &operator*=(const FTps<T, N> &rhs);
    
    /// Multiply.
    FVps operator*(const FVps<T, N>& rhs) const;

    /// Divide and assign.
    //  Throw DivideError if constant part of [b]rhs[/b] is zero.
    FVps &operator/=(const FTps<T, N> &rhs);

    /// Multiply and assign.
    FVps &operator*=(const T &rhs);

    /// Divide and assign.
    //  Throw DivideError if [b]rhs[/b] is zero.
    FVps &operator/=(const T &rhs);

    /// Inverse.
    //  Invert a map.
    FVps inverse(int trunc = (FTps<T, N>::EXACT)) const;

    /// Inverse.
    //  Invert a map.
    FVps myInverse(int trunc = (FTps<T, N>::EXACT)) const;

    /// Partial derivative.
    //  Return partial derivative with respect to variable [b]var[/b].
    FVps derivative(int var) const;

    /// Partial integral.
    //  Return partial integral with respect to variable [b]var[/b].
    //  Throw LogicalError for a constant.
    FVps integral(int var) const;

    /// Extract the constant part of the map.
    //  This is equivalent to evaluating the map at the origin.
    FVector<T, N> constantTerm() const;

    /// Extract the constant part of the map expanded about point [b]P[/b].
    //  This is equivalent to evaluating the map at point [b]P[/b].
    FVector<T, N> constantTerm(const FVector<T, N> &P) const;

    /// Extract the linear part of the map.
    FMatrix<T, N, N> linearTerms() const;

    /// Extract the linear part of the map expanded about point [b]P[/b].
    //  This is equivalent to the Jacobian matrix at point [b]P[/b].
    FMatrix<T, N, N> linearTerms(const FVector<T, N> &P) const;

    /// Return orders {min, max, trc} of f(rhs(z)).
    // NB: This routine is NOT guaranteed to return min <= max <= trc.
    // If min exceeds max, then f(rhs(z)) = 0 + O(z^{trc+1}).
    // This possibility must be checked for.
    Array1D<int> getSubstOrders(const FVps<T, N> &rhs, int trunc = (FTps<T, N>::EXACT)) const;

    /// Substitute.
    //  Use the linear map with matrix representation [b]M[/b] to transform the
    //  order n part of an FVps; leave remaining parts unchanged.  Return a new FVps.
    // NB: This method uses static local memory.
    FVps substitute(const FMatrix<T, N, N> &M, int n) const;

    /// Substitute.
    //  Use the linear map with matrix representation [b]M[/b] to transform the
    //  order nl through nh parts of an FVps; leave remaining parts unchanged.
    //  Return a new FVps.
    // NB: This method uses static local memory.
    FVps substitute(const FMatrix<T, N, N> &M, int nl, int nh) const;

    /// Substitute.
    //  Use the linear map with matrix representation [b]M[/b]
    //  to transform an FVps.  Return a new FVps.
    // NB: This method (indirectly) uses static local memory.
    FVps substitute(const FMatrix<T, N, N> &M) const;

    /// Substitute.
    //  Use the (nonlinear) map [b]m[/b] to transform an FVps.  Return a new FVps.
    // NB: This method uses static local memory.
    FVps substitute(const FVps<T, N> &m, int trunc = (FTps<T, N>::EXACT)) const;

    /// Substitute map into matrix.
    FVps substituteInto(const FMatrix<T, N, N> &lhs) const;
    
    /// Get a FTps that is a combination of the polynomials of FVps.
    // Computes a FTps by multiplying the FTps of FVps using the powers specified by [b]power[/b].
    FTps<T, N> getFTps(const FArray1D<int, N>& power) const;
    
    /// Get a FVps from stream [b]is[/b].
    std::istream &get(std::istream &is);

    /// Put a FVps to stream [b]os[/b].
    std::ostream &put(std::ostream &os) const;

protected:

    // Representation of the FVps.
    FTps<T, N> data[N];

};


// Global functions.
// ------------------------------------------------------------------------

/// Add.
template <class T, int N>
FVps<T, N> operator+(const FVps<T, N> &lhs, const FVps<T, N> &rhs);

/// Subtract.
template <class T, int N>
FVps<T, N> operator-(const FVps<T, N> &lhs, const FVps<T, N> &rhs);

/// Add.
template <class T, int N>
FVps<T, N> operator+(const FVps<T, N> &lhs, const FVector<T, N> &rhs);

/// Subtract.
template <class T, int N>
FVps<T, N> operator-(const FVps<T, N> &lhs, const FVector<T, N> &rhs);

/// Add.
template <class T, int N>
FVps<T, N> operator+(const FVector<T, N> &lhs, const FVps<T, N> &rhs);

/// Subtract.
template <class T, int N>
FVps<T, N> operator-(const FVector<T, N> &lhs, const FVps<T, N> &rhs);

/// Multiply.
template <class T, int N>
FVector<T, N> operator*(const FVps<T, N> &lhs, const FVector<T, N>& rhs);

/// Multiply.
template <class T, int N>
FVps<T, N> operator*(const FVps<T, N> &lhs, const FTps<T, N> &rhs);

/// Multiply.
template <class T, int N>
FVps<T, N> operator*(const FTps<T, N> &lhs, const FVps<T, N> &rhs);

/// Multiply.
template <class T, int N>
FVps<T, N> operator*(const FVps<T, N> &lhs, const T &rhs);

/// Multiply.
template <class T, int N>
FVps<T, N> operator*(const T &lhs, const FVps<T, N> &rhs);

/// Multiply.
template <class T, int N>
FVps<T, N> operator*(const FMatrix<T, N, N> &lhs, const FVps<T, N> &rhs);

/// Divide.
//  Throw DivideError, if constant part of [b]rhs[/b] is zero.
template <class T, int N>
FVps<T, N> operator/(const FVps<T, N> &lhs, const FTps<T, N> &rhs);

/// Divide.
//  Throw DivideError, if [b]rhs[/b] is zero.
template <class T, int N>
FVps<T, N> operator/(const FVps<T, N> &lhs, const T &rhs);

/// Build the exponential series.
//  Return the series exp(:H:) M,
//  the Lie transform exp(:H:) acting on the map [b]M[/b].
template <class T, int N> FVps<T, N>
ExpMap(const FTps<T, N> &H, const FVps<T, N> &M, int trunc = (FTps<T, N>::EXACT));

/// Poisson bracket.
template <class T, int N> FVps<T, N>
PoissonBracket(const FTps<T, N> &x, const FVps<T, N> &y, int trunc = (FTps<T, N>::EXACT));

/// Extract FVps from stream [b]is[/b].
template <class T, int N>
std::istream &operator>>(std::istream &is, FVps<T, N> &);

/// Insert FVps to stream [b]os[/b].
template <class T, int N>
std::ostream &operator<<(std::ostream &os, const FVps<T, N> &vps);


// Implementation.
#include "FixedAlgebra/FVps.hpp"

#endif // CLASSIC_FVps_HH
