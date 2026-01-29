#ifndef CLASSIC_Vps_HH
#define CLASSIC_Vps_HH

// ------------------------------------------------------------------------
// $RCSfile: Vps.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Vps
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

template <class T> class Matrix;
template <class T> class Vps;

#include "Algebra/Tps.h"
#include "Algebra/Vector.h"
#include <iosfwd>


// Template class Vps<T>
// ------------------------------------------------------------------------
/// Vector truncated power series.
//  A vector of truncated power series with coefficients of type [b]T[/b]
//  in [b]n[/b] variables. This class serves as the base class for maps.

template <class T>
class Vps {

public:

    /// Constructor.
    //  Construct [b]nDim[/b] series in [b]nVar[/b] variables each.
    Vps(int nDim, int nVar = 0);

    /// Convert.
    //  The constant part is zero.
    //  The linear part is set to the contents of [b]M[/b].
    Vps(const Matrix<T> &M);

    /// Convert from vector.
    //  The constant part is set to the contents of [b]V[/b].
    //  The linear part is zero.
    Vps(const Vector<T> &V);

    Vps();
    Vps(const Vps<T> &rhs);
    ~Vps();
    Vps<T> &operator=(const Vps<T> &);

    /// Get component.
    //  Return the component [b]index[/b].
    //  Throw RangeError, if [b]index[/b] is out of range.
    const Tps<T> &getComponent(int index) const;

    /// Set component.
    //  Set the component [b]index[/b].
    //  Throw RangeError, if [b]index[/b] is out of range.
    void setComponent(int index, const Tps<T> &value);

    /// Get component.
    //  Return the component [b]index[/b].
    //  Result is undefined, if [b]index[/b] is out of range.
    Tps<T> &operator[](int index);

    /// Set component.
    //  Set the component [b]index[/b].
    //  Result is undefined, if [b]index[/b] is out of range.
    const Tps<T> &operator[](int index) const;

    /// Check consistency.
    void check() const;

    /// Unary plus.
    Vps<T> operator+() const;

    /// Unary minus.
    Vps<T> operator-() const;

    /// Multiply by Tps<T> and assign.
    Vps<T> &operator*=(const Tps<T> &y);

    /// Divide by Tps<T> and assign.
    // Throw DivideError if constant part of [b]y[/b] is zero.
    Vps<T> &operator/=(const Tps<T> &y);

    /// Multiply by constant and assign.
    Vps<T> &operator*=(const T &y);

    /// Divide by constant and assign.
    //  Throw DivideError if [b]y[/b] is zero.
    Vps<T> &operator/=(const T &y);

    /// Addition.
    Vps<T> &operator+=(const Vps<T> &y);

    /// Subtraction.
    Vps<T> &operator-=(const Vps<T> &y);

    /// Add and assign.
    Vps<T> &operator+=(const Vector<T> &y);

    /// Subtract and assign.
    Vps<T> &operator-=(const Vector<T> &y);

    /// Get a Vps<T> from stream is.
    std::istream &get(std::istream &is);

    /// Put a Vps<T> to stream os.
    std::ostream &put(std::ostream &os) const;

    /// Get dimension (number of Tps<T> components).
    int getDimension() const;

    /// Get highest order contained in any component.
    int getTopOrder() const;

    /// Get lowest truncation order in any component.
    int getTruncOrder() const;

    /// Get number of variables (the same in all components).
    int getVariables() const;

    /// Extract range of orders, set others to zero.
    Vps<T> filter(int lowOrder, int highOrder) const;

    /// Truncate, may also increase truncation order.
    Vps<T> truncate(int trunc);

protected:

    // Representation of the Vps<T>.
    Array1D< Tps<T> > data;

    // The number of variables.
    mutable int variables;
};


// Unchecked inline access to components.
// ------------------------------------------------------------------------

template <class T> inline
Tps<T> &Vps<T>::operator[](int index)
{ return data[index]; }

template <class T> inline
const Tps<T> &Vps<T>::operator[](int index) const
{ return data[index]; }


// Global Operators on Vps<T>
// ------------------------------------------------------------------------

/// Multiply.
template <class T> Vps<T> operator*(const Vps<T> &x, const Tps<T> &y);

/// Multiply.
template <class T> Vps<T> operator*(const Tps<T> &x, const Vps<T> &y);

/// Multiply.
template <class T> Vps<T> operator*(const Vps<T> &x, const T &y);

/// Multiply.
template <class T> Vps<T> operator*(const T &x, const Vps<T> &y);

/// Divide.
//  Throw DivideError, if constant part of [b]y[/b] is zero.
template <class T> Vps<T> operator/(const Vps<T> &x, const Tps<T> &y);

/// Divide.
//  Throw DivideError, if [b]y[/b] is zero.
template <class T> Vps<T> operator/(const Vps<T> &x, const T &y);

/// Add.
template <class T> Vps<T> operator+(const Vps<T> &x, const Vps<T> &y);

/// Subtract.
template <class T> Vps<T> operator-(const Vps<T> &x, const Vps<T> &y);

/// Add.
template <class T> Vps<T> operator+(const Vps<T> &x, const Vector<T> &y);

/// Subtract.
template <class T> Vps<T> operator-(const Vps<T> &x, const Vector<T> &y);

/// Add.
template <class T> Vps<T> operator+(const Vector<T> &x, const Vps<T> &y);

/// Subtract.
template <class T> Vps<T> operator-(const Vector<T> &x, const Vps<T> &y);

/// Extract from stream.
template <class T> std::istream &operator>>(std::istream &, Vps<T> &x);

/// Insert to stream.
template <class T> std::ostream &operator<<(std::ostream &, const Vps<T> &x);


// Implementation of global functions for class Vps<T>.
// ------------------------------------------------------------------------

template <class T> inline
Vps<T> operator*(const Vps<T> &x, const Tps<T> &y)
{ Vps<T> z(x); return z *= y; }

template <class T> inline
Vps<T> operator*(const Tps<T> &x, const Vps<T> &y)
{ Vps<T> z(x); return z *= y; }

template <class T> inline
Vps<T> operator*(const Vps<T> &x, const T &y)
{ Vps<T> z(x); return z *= y; }

template <class T> inline
Vps<T> operator*(const T &x, const Vps<T> &y)
{ Vps<T> z(x); return z *= y; }

template <class T> inline
Vps<T> operator/(const Vps<T> &x, const Tps<T> &y)
{ Vps<T> z(x); return z /= y; }

template <class T> inline
Vps<T> operator/(const Vps<T> &x, const T &y)
{ Vps<T> z(x); return z /= y; }


template <class T> inline
Vps<T> operator+(const Vps<T> &x, const Vps<T> &y)
{ Vps<T> z(x); return z += y; }

template <class T> inline
Vps<T> operator-(const Vps<T> &x, const Vps<T> &y)
{ Vps<T> z(x); return z -= y; }

template <class T> inline
Vps<T> operator+(const Vps<T> &x, const Vector<T> &y)
{ Vps<T> z(x); return z += y; }

template <class T> inline
Vps<T> operator-(const Vps<T> &x, const Vector<T> &y)
{ Vps<T> z(x); return z -= y; }

template <class T> inline
Vps<T> operator+(const Vector<T> &x, const Vps<T> &y)
{ Vps<T> z(y); return z += x; }

template <class T> inline
Vps<T> operator-(const Vector<T> &x, const Vps<T> &y)
{ Vps<T> z(- y); return z += x; }

template <class T> inline
std::istream &operator>>(std::istream &is, Vps<T> &x)
{ return x.get(is); }

template <class T> inline
std::ostream &operator<<(std::ostream &os, const Vps<T> &x)
{ return x.put(os); }

#endif // CLASSIC_Vps_HH
