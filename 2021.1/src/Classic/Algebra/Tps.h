#ifndef CLASSIC_Tps_HH
#define CLASSIC_Tps_HH

// ------------------------------------------------------------------------
// $RCSfile: Tps.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: Tps<class T>
//   Truncated power series in n variables of type T.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2002/03/25 20:44:14 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include <iosfwd>
#include <climits>

template <class T> class Tps;
template <class T> class TpsRep;
template <class T> class Matrix;
template <class T> class Vector;
template <class T> class VpsMap;

class TpsMonomial;


// Class Tps
// ------------------------------------------------------------------------
/// Truncated power series.
//  Truncated power series in [b]n[/b] variables of type [b]T[/b].
//  All divide operations throw DivideError,
//  if the constant part of the divisor is zero.
//  All operations throw SizeError, if the operands are not consistent.
//  All operations throw RangeError, if an index is out of range.

template <class T>
class Tps {

public:

    /// Constructor.
    //  Define maximum order and number of variables, value = 0.
    Tps(int maxOrder, int nVar);

    /// Conversion.
    //  From constant [b]y[/b].
    Tps(const T &y);

    /// Conversion.
    //  From int [b]y[/b].
    Tps(int y);

    Tps();
    Tps(const Tps<T> &x);
    ~Tps();
    Tps<T> &operator=(const Tps<T> &y);

    /// Convert and assign.
    Tps<T> &operator=(const T &y);

    /// Extract orders.
    //  Return the orders [b]lowOrder[/b] through [b]highOrder[/b],
    //  both included.
    Tps<T> filter(int lowOrder, int highOrder) const;

    /// Truncate.
    //  Change the maximum order to [b]trunc[/b]; may also reserve more space.
    Tps<T> truncate(int trunc);

    /// Get coefficient.
    //  Return the coefficient denoted by the Giorgilli index.
    //  Throw RangeError, if index is out of range.
    const T getCoefficient(int index) const;

    /// Set coefficient.
    //  Assign the coefficient denoted by the Giorgilli index.
    //  Expand size if index is out of range.
    void setCoefficient(int index, const T &value);

    /// Get coefficient.
    //  Return the coefficient denoted by monomial exponents.
    //  Throw RangeError, if index is out of range.
    const T getCoefficient(const TpsMonomial &monomial) const;

    /// Set coefficient.
    //  Assign the coefficient denoted by monomial exponents.
    //  Throw RangeError, if index is out of range.
    void setCoefficient(const TpsMonomial &monomial, const T &value);

    /// Get coefficient.
    //  Return the coefficient denoted by the Giorgilli index.
    //  Result is undefined, if index is out of range.
    const T operator[](int index) const;

    /// Get coefficient.
    //  REturn a reference to the coefficient denoted by the Giorgilli index.
    //  Result is undefined, if index is out of range.
    T &operator[](int index);

    /// Get coefficient.
    //  Return the coefficient denoted by monomial exponents.
    //  Result is undefined, if index is out of range.
    const T operator[](const TpsMonomial &monomial) const;

    /// Get coefficient.
    //  Return reference to the coefficient denoted by monomial exponents.
    //  Result is undefined, if index is out of range.
    T &operator[](const TpsMonomial &monomial);

    /// Make variable.
    //  Construct the variable identified by the index [b]var[/b],
    //  with total of [b]nVar[/b] variables.
    static Tps<T> makeVariable(int nVar, int var);

    /// Make power.
    //  Construct power of variable [b]var[/b], with total of [b]nVar[/b]
    //  variables.
    static Tps<T> makeVarPower(int nVar, int var, int power);

    /// Make monomial.
    //  Construct the monomial with the exponents in [b]m[/b] and coefficient
    //  [b]t[/b].
    static Tps<T> makeMonomial(const TpsMonomial &m, const T &t);

    /// Unary plus.
    Tps<T> operator+() const;

    /// Unary minus.
    Tps<T> operator-() const;

    /// Add and assign.
    Tps<T> &operator+=(const Tps<T> &y);

    /// Subtract and assign.
    Tps<T> &operator-=(const Tps<T> &y);

    /// Multiply and assign.
    Tps<T> &operator*=(const Tps<T> &y);

    /// Divide and assign.
    Tps<T> &operator/=(const Tps<T> &y);

    /// Add constant and assign.
    Tps<T> &operator+=(const T &y);

    /// Subtract constant and assign.
    Tps<T> &operator-=(const T &y);

    /// Multiply by constant and assign.
    Tps<T> &operator*=(const T &y);

    /// Divide by constant and assign.
    Tps<T> &operator/=(const T &y);

    /// Equality operator.
    bool operator==(const Tps<T> &y) const;

    /// Equlatity with constant.
    bool operator==(const T &y) const;

    /// Inequality operator.
    bool operator!=(const Tps<T> &y) const;

    /// Inequality with constant.
    bool operator!=(const T &y) const;

    /// Substitute.
    //  Substitute matrix [b]M[/b] in Tps, giving new Tps.
    Tps<T> substitute(const Matrix<T> &M) const;

    /// Substitute.
    //  Substitute map [b]m[/b] in Tps, giving new Tps.
    Tps<T> substitute(const VpsMap<T> &m) const;

    /// Substitute.
    //  Substitute vector [b]v[/b] in Tps, giving a constant.
    T evaluate(const Vector<T> &v) const;

    /// Set to zero.
    void clear();

    /// Reciprocal value.
    Tps<T> inverse(int order = truncOrder) const;

    /// Truncated multiplication.
    Tps<T> multiply(const Tps<T> &y, int trunc) const;

    /// Partial derivative.
    //  Return partial derivative with respect to variable [b]var[/b].
    //  Return zero for a constant.
    Tps<T> derivative(int var) const;

    /// Partial integral.
    //  Return partial integral with respect to variable [b]var[/b].
    //  Throw LogicalError for a constant.
    Tps<T> integral(int var) const;

    /// Multiply by variable [b]var[/b].
    //  Throw LogicalError for a constant.
    Tps<T> multiplyVariable(int var) const;

    /// Multiply monomial-wise.
    //  Construct new Tps<T> by multipying corresponding monomials pair-wise.
    Tps<T> scaleMonomials(const Tps<T> &y) const;

    /// Taylor series.
    //  Expand truncated Taylor series with coefficiens [b]series[/b]
    //  and order [b]n[/b].
    Tps<T> Taylor(const T series[], int n) const;

    /// Get maximal order.
    int getMaxOrder() const;

    /// Get truncation order.
    int getTruncOrder() const;

    /// Get number of variables.
    int getVariables() const;

    /// Get number of coefficients.
    int getSize() const;

    /// Get Tps from the stream is.
    std::istream &get(std::istream &is);

    /// Put Tps to the stream is.
    std::ostream &put(std::ostream &os) const;

    /// Get exponents.
    //  Return the exponent array for the given Giorgilli index.
    const TpsMonomial &getExponents(int index) const;

    /// Get order.
    //  Return the order of the monomial with the given Giorgilli index.
    int getOrder(int index) const;

    /// Size.
    //  Return the number of monomials required for a Tps of order [b]order[/b].
    int getSize(int order) const;

    /// Test for constant.
    bool isConstant() const;

    /// Get global truncation order.
    static int getGlobalTruncOrder();

    /// Set global truncation order.
    static void setGlobalTruncOrder(int order);

    /// Representation of infinite precision.
    static const int EXACT = INT_MAX;

private:

    // Internal constructor.
    Tps(TpsRep<T> *);

    // Set maximum order.
    void setMaxOrder(int);

    // Make representation unique.
    void unique();

    // Pointer to representation.
    TpsRep<T> *rep;

    // Current maximum truncation order.
    static int truncOrder;
};


// Global functions acting on Tps<T> objects.
// ------------------------------------------------------------------------

/// Add.
template <class T> Tps<T> operator+(const Tps<T> &x, const Tps<T> &y);

/// Subtract.
template <class T> Tps<T> operator-(const Tps<T> &x, const Tps<T> &y);

/// Add.
template <class T> Tps<T> operator+(const Tps<T> &x, const T &y);

/// Subtract.
template <class T> Tps<T> operator-(const Tps<T> &x, const T &y);

/// Add.
template <class T> Tps<T> operator+(const T &x, const Tps<T> &y);

/// Subtract.
template <class T> Tps<T> operator-(const T &x, const Tps<T> &y);

/// Multiply.
template <class T> Tps<T> operator*(const Tps<T> &x, const Tps<T> &y);

/// Divide.
template <class T> Tps<T> operator/(const Tps<T> &x, const Tps<T> &y);

/// Multiply.
template <class T> Tps<T> operator*(const Tps<T> &x, const T &y);

/// Divide.
template <class T> Tps<T> operator/(const Tps<T> &x, const T &y);

/// Multiply.
template <class T> Tps<T> operator*(const T &x, const Tps<T> &y);

/// Divide.
template <class T> Tps<T> operator/(const T &x, const Tps<T> &y);

/// Test for equality.
template <class T> bool operator==(const T &x, const Tps<T> &y);

/// Extract from stream.
template <class T> std::istream &operator>>(std::istream &, Tps<T> &x);

/// Insert to stream.
template <class T> std::ostream &operator<<(std::ostream &, const Tps<T> &x);


// Inline global functions.
// ------------------------------------------------------------------------

template <class T> Tps<T> operator+(const Tps<T> &x, const Tps<T> &y)
{ Tps<T> z(x); return z += y; }

template <class T> Tps<T> operator-(const Tps<T> &x, const Tps<T> &y)
{ Tps<T> z(x); return z -= y; }

template <class T> inline
Tps<T> operator+(const Tps<T> &x, const T &y)
{ Tps<T> z(x); return z += y; }

template <class T> inline
Tps<T> operator-(const Tps<T> &x, const T &y)
{ Tps<T> z(x); return z -= y; }

template <class T> inline
Tps<T> operator+(const T &x, const Tps<T> &y)
{ Tps<T> z(y); return z += x; }

template <class T>
Tps<T> operator-(const T &x, const Tps<T> &y)
{ Tps<T> z(- y); return z += x; }

template <class T> Tps<T> operator*(const Tps<T> &x, const Tps<T> &y)
{ Tps<T> z(x); return z *= y; }

template <class T> Tps<T> operator/(const Tps<T> &x, const Tps<T> &y)
{ Tps<T> z(x); return z /= y; }

template <class T> Tps<T> operator*(const Tps<T> &x, const T &y)
{ Tps<T> z(x); return z *= y; }

template <class T> Tps<T> operator/(const Tps<T> &x, const T &y)
{ Tps<T> z(x); return z /= y; }

template <class T> Tps<T> operator*(const T &x, const Tps<T> &y)
{ Tps<T> z(y); return z *= x; }

template <class T> Tps<T> operator/(const T &x, const Tps<T> &y)
{ Tps<T> z(y.inverse(Tps<T>::getGlobalTruncOrder())); return z *= x; }

template <class T> bool operator==(const T &x, const Tps<T> &y)
{ return y == x; }

template <class T> std::istream &operator>>(std::istream &is, Tps<T> &x)
{ return x.get(is); }

template <class T> std::ostream &operator<<(std::ostream &os, const Tps<T> &x)
{ return x.put(os); }

#endif // CLASSIC_Tps_HH
