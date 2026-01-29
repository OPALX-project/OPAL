#ifndef CLASSIC_LinearFun_HH
#define CLASSIC_LinearFun_HH

// ------------------------------------------------------------------------
// $RCSfile: LinearFun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: LinearFun<T,N>
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
template <class T, int N> class LinearMap;
template <class T, int M, int N> class FMatrix;


// Template class LinearFun<T,N>.
// ------------------------------------------------------------------------
/// Linear function in N variables of type T.
//  All divide operations throw DivideError,
//  if the constant part of the divisor is zero.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler do the correct thing, and are not defined for speed.

template <class T, int N>
class LinearFun {

public:

    /// Default constructor.
    //  Construct zero value.
    LinearFun();

    /// Conversion.
    LinearFun(const T &);

    /// Conversion.
    LinearFun(int);

    /// Convert and assign.
    LinearFun &operator=(const T &y);

    /// Get coefficient.
    //  Return value of the coefficient denoted by the index.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Return zero, if index is out of range.
    const T getCoefficient(int index) const;

    /// Set coefficient.
    //  Assign value of the coefficient denoted by the index.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Ignore, if index is out of range.
    void setCoefficient(int index, const T &value);

    /// Get coefficient.
    //  Return value of the coefficient denoted by the index.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Result is undefined, if index is out of range.
    inline const T operator[](int index) const;

    /// Get coefficient.
    //  Return a reference to the coefficient denoted by the index.
    //  Result is undefined, if index is out of range.
    //  The index zero stands for the constant term,
    //  any other value selects the derivative by the variable index+1.
    //  Use as a lvalue is only allowed when the LinearFun is unique().
    inline T &operator[](int index);

    /// Make variable.
    //  Construct the variable identified by the index [b]var[/b],
    //  with total of [b]nVar[/b] variables.
    static LinearFun makeVariable(int var);

    /// Unary plus.
    LinearFun operator+() const;

    /// Unary minus.
    LinearFun operator-() const;

    /// Add and assign.
    LinearFun &operator+=(const LinearFun &y);

    /// Subtract and assign.
    LinearFun &operator-=(const LinearFun &y);

    /// Multiply and assign.
    LinearFun &operator*=(const LinearFun &y);

    /// Approximate division and assignation.
    LinearFun &operator/=(const LinearFun &y);

    /// Add constant and assign.
    LinearFun &operator+=(const T &y);

    /// Subtract constant and assign.
    LinearFun &operator-=(const T &y);

    /// Multiply by constant and assign.
    LinearFun &operator*=(const T &y);

    /// Divide by constant and assign.
    LinearFun &operator/=(const T &y);

    /// Equality operator.
    bool operator==(const LinearFun &y) const;

    /// Equality with constant.
    bool operator==(const T &y) const;

    /// Inequality operator.
    bool operator!=(const LinearFun &y) const;

    /// Inequality with constant.
    bool operator!=(const T &y) const;

    /// Approximate reciprocal value 1/(*this).
    LinearFun inverse() const;

    /// Multiplication truncated to order one.
    LinearFun multiply(const LinearFun &y) const;

    /// Evaluate LinearFun at point.
    T evaluate(const FVector<T, N> &) const;

    /// Substitute.
    //  Substitute map [b]m[/b] in LinearFun, giving new LinearFun.
    LinearFun<T, N> substitute(const LinearMap<T, N> &m) const;

    /// Substitute.
    //  Substitute matrix [b]M[/b] in LinearFun, giving new LinearFun.
    LinearFun<T, N> substitute(const FMatrix<T, N, N> &M) const;

    /// Taylor series.
    //  Expand Taylor series with coefficiens [b]series[/b] and order one.
    LinearFun taylor(const T series[2]) const;

    /// Read LinearFun on the stream [b]is[/b].
    std::istream &get(std::istream &is);

    /// Write LinearFun on the stream [b]os[/b].
    std::ostream &put(std::ostream &os) const;

private:

    /// Representation.
    //  data[0] is the constant term,
    //  data[v+1] is the coefficient of variable [b]v[/b].
    T data[N+1];
};


// Global functions.
// ------------------------------------------------------------------------

/// Add.
template <class T, int N>
LinearFun<T, N> operator+(const LinearFun<T, N> &, const LinearFun<T, N> &);

/// Subtract.
template <class T, int N>
LinearFun<T, N> operator-(const LinearFun<T, N> &, const LinearFun<T, N> &);

/// Add.
template <class T, int N>
LinearFun<T, N> operator+(const LinearFun<T, N> &, const T &);

/// Subtract.
template <class T, int N>
LinearFun<T, N> operator-(const LinearFun<T, N> &, const T &);

/// Add.
template <class T, int N>
LinearFun<T, N> operator+(const T &, const LinearFun<T, N> &);

/// Subtract.
template <class T, int N>
LinearFun<T, N> operator-(const T &, const LinearFun<T, N> &);

/// Multiply.
template <class T, int N>
LinearFun<T, N> operator*(const LinearFun<T, N> &, const LinearFun<T, N> &);

/// Divide.
template <class T, int N>
LinearFun<T, N> operator/(const LinearFun<T, N> &, const LinearFun<T, N> &);

/// Multiply.
template <class T, int N>
LinearFun<T, N> operator*(const LinearFun<T, N> &, const T &);

/// Divide.
template <class T, int N>
LinearFun<T, N> operator/(const LinearFun<T, N> &, const T &);

/// Multiply.
template <class T, int N>
LinearFun<T, N> operator*(const T &, const LinearFun<T, N> &);

/// Divide 1 / LinearFun.
template <class T, int N>
LinearFun<T, N> operator/(const T &, const LinearFun<T, N> &);

/// Equality.
template <class T, int N>
bool operator==(const T &, const LinearFun<T, N> &);

/// Inequality.
template <class T, int N>
bool operator!=(const T &, const LinearFun<T, N> &);

/// Extract LinearFun from stream [b]is[/b].
template <class T, int N>
std::istream &operator>>(std::istream &is, LinearFun<T, N> &);

/// Insert LinearFun to stream [b]os[/b].
template <class T, int N>
std::ostream &operator<<(std::ostream &os, const LinearFun<T, N> &);


// Implementation.
#include "FixedAlgebra/LinearFun.hpp"
#include "FixedAlgebra/LinearMap.h"

#endif // CLASSIC_LinearFun_HH
