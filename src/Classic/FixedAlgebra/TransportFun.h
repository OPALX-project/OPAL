#ifndef CLASSIC_TransportFun_HH
#define CLASSIC_TransportFun_HH

// ------------------------------------------------------------------------
// $RCSfile: TransportFun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: TransportFun<T,N>
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

template <class T, int N> class FVector;
template <class T, int N> class TransportMap;
template <class T, int M, int N> class FMatrix;


// Template class TransportFun<T,N>.
// ------------------------------------------------------------------------
/// Transport function in N variables of type T.
//  All divide operations throw DivideError,
//  if the constant part of the divisor is zero.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons
//  they are not implemented.

template <class T, int N>
class TransportFun {

public:

    /// Default constructor.
    //  Construct zero value.
    TransportFun();

    /// Conversion.
    TransportFun(const T &);

    /// Conversion.
    TransportFun(int);

    /// Convert and assign.
    TransportFun &operator=(const T &y);

    /// Get coefficient.
    //  Return value of the coefficient denoted by the index.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Return zero, if index is out of range.
    const T getCoefficient(int index) const;

    /// Get coefficient.
    //  Return value of the second-order coefficient denoted by indices.
    //  Return zero, if indices are out of range.
    const T getCoefficient(int index1, int index2) const;

    /// Set coefficient.
    //  Assign value of the coefficient denoted by the index.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Ignore, if index is out of range.
    void setCoefficient(int index, const T &value);

    /// Set coefficient.
    //  Assign value of the second-order coefficient denoted by indices.
    //  Ignore, if indices are out of range.
    void setCoefficient(int index1, int index2, const T &value);

    /// Get coefficient.
    //  Return value of the coefficient denoted by the index.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Result is undefined, if index is out of range.
    inline const T operator[](int) const;

    /// Get coefficient.
    //  Return a reference to the coefficient denoted by the index.
    //  Result is undefined, if index is out of range.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    inline T &operator[](int);

    /// Get coefficient.
    //  Return value of the second-order coefficient denoted by indices.
    //  NOTE: i1 <= i2.
    //  Result is undefined, if indices are out of range.
    inline const T operator()(int i1, int i2) const;

    /// Get coefficient.
    //  Return a reference to the second-order coefficient denoted by indices.
    //  NOTE: i1 <= i2.
    //  Result is undefined, if indices are out of range.
    inline T &operator()(int i1, int i2);

    /// Make variable.
    //  Construct the variable identified by the index [b]var[/b],
    //  with total of [b]nVar[/b] variables.
    static TransportFun makeVariable(int var);

    /// Unary plus.
    TransportFun operator+() const;

    /// Unary minus.
    TransportFun operator-() const;

    /// Add and assign.
    TransportFun &operator+=(const TransportFun &y);

    /// Subtract and assign.
    TransportFun &operator-=(const TransportFun &y);

    /// Multiply and assign.
    TransportFun &operator*=(const TransportFun &y);

    /// Approximate division and assignation.
    TransportFun &operator/=(const TransportFun &y);

    /// Add constant and assign.
    TransportFun &operator+=(const T &y);

    /// Subtract constant and assign.
    TransportFun &operator-=(const T &y);

    /// Multiply by constant and assign.
    TransportFun &operator*=(const T &y);

    /// Divide by constant and assign.
    TransportFun &operator/=(const T &y);

    /// Equality operator.
    bool operator==(const TransportFun &y) const;

    /// Equality with constant.
    bool operator==(const T &y) const;

    /// Inequality operator.
    bool operator!=(const TransportFun &y) const;

    /// Inequality with constant.
    bool operator!=(const T &y) const;

    /// Approximate reciprocal value 1/(*this).
    TransportFun inverse() const;

    /// Multiplication truncated to order one.
    TransportFun multiply(const TransportFun &y) const;

    /// Multiply by variable var.
    TransportFun<T, N> multiplyVariable(int var) const;

    /// Evaluate TransportFun at point.
    T evaluate(const FVector<T, N> &) const;

    /// Substitute.
    //  Substitute map [b]m[/b] in TransportFun, giving new TransportFun.
    TransportFun<T, N> substitute(const TransportMap<T, N> &m) const;

    /// Substitute.
    //  Substitute matrix [b]M[/b] in TransportFun, giving new TransportFun.
    TransportFun<T, N> substitute(const FMatrix<T, N, N> &M) const;

    /// Taylor series.
    //  Expand Taylor series with coefficiens [b]series[/b] and order one.
    TransportFun taylor(const T series[2]) const;

    /// Read TransportFun on the stream [b]is[/b].
    std::istream &get(std::istream &is);

    /// Write TransportFun on the stream [b]os[/b].
    std::ostream &put(std::ostream &os) const;

private:

    /// Size of representation.
    static const int SIZE = (N + 1) *(N + 2) / 2;

    /// Representation.
    //  data[0] is the constant term,
    //  data[v+1] is the coefficient of variable [b]v[/b].
    T data[SIZE];
};


// Global functions.
// ------------------------------------------------------------------------

/// Add.
template <class T, int N>
TransportFun<T, N> operator+(const TransportFun<T, N> &, const TransportFun<T, N> &);

/// Subtract.
template <class T, int N>
TransportFun<T, N> operator-(const TransportFun<T, N> &, const TransportFun<T, N> &);

/// Add.
template <class T, int N>
TransportFun<T, N> operator+(const TransportFun<T, N> &, const T &);

/// Subtract.
template <class T, int N>
TransportFun<T, N> operator-(const TransportFun<T, N> &, const T &);

/// Add.
template <class T, int N>
TransportFun<T, N> operator+(const T &, const TransportFun<T, N> &);

/// Subtract.
template <class T, int N>
TransportFun<T, N> operator-(const T &, const TransportFun<T, N> &);

/// Multiply.
template <class T, int N>
TransportFun<T, N> operator*(const TransportFun<T, N> &, const TransportFun<T, N> &);

/// Divide.
template <class T, int N>
TransportFun<T, N> operator/(const TransportFun<T, N> &, const TransportFun<T, N> &);

/// Multiply.
template <class T, int N>
TransportFun<T, N> operator*(const TransportFun<T, N> &, const T &);

/// Divide.
template <class T, int N>
TransportFun<T, N> operator/(const TransportFun<T, N> &, const T &);

/// Multiply.
template <class T, int N>
TransportFun<T, N> operator*(const T &, const TransportFun<T, N> &);

/// Divide 1 / TransportFun.
template <class T, int N>
TransportFun<T, N> operator/(const T &, const TransportFun<T, N> &);

/// Equality.
template <class T, int N>
bool operator==(const T &, const TransportFun<T, N> &);

/// Inequality.
template <class T, int N>
bool operator!=(const T &, const TransportFun<T, N> &);

/// Extract TransportFun from stream [b]is[/b].
template <class T, int N>
std::istream &operator>>(std::istream &is, TransportFun<T, N> &);

/// Insert TransportFun to stream [b]os[/b].
template <class T, int N>
std::ostream &operator<<(std::ostream &os, const TransportFun<T, N> &);


// Implementation.
#include "FixedAlgebra/TransportFun.hpp"
#include "FixedAlgebra/TransportMap.h"

#endif // CLASSIC_TransportFun_HH
