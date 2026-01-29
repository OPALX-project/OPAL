#ifndef CLASSIC_Taylor_HH
#define CLASSIC_Taylor_HH

// ------------------------------------------------------------------------
// $RCSfile: Taylor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Taylor<T>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2001/11/15 08:52:26 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"


// Template class Taylor<T>.
// ------------------------------------------------------------------------
/// A representation for a Taylor series in one variable,
//  whose coefficients are of type T.  For some operations this class
//  assumes that the type T has particular operations, like PoissonBracket.
//  The class Taylor<T> is mainly meant for implementing a Lie generator
//  with one parameter as a Taylor series in this parameter.  This permits
//  easy integration with respect to the parameter.

template <class T>
class Taylor {

public:

    /// Construct a zero series of a given order.
    explicit Taylor(int);

    Taylor();
    Taylor(const Taylor &);
    ~Taylor();
    const Taylor &operator=(const Taylor &);


    /// Get pointer to beginning of series (zero-order term).
    //  Version for non-constant series.
    inline T *begin();

    /// Get pointer to beginning of series (zero-order term).
    //  Version for constant series.
    inline const T *begin() const;

    /// Get pointer to end of series (one beyond highest term).
    //  Version for non-constant series.
    inline T *end();

    /// Get pointer to end of series (one beyond highest term).
    //  Version for constant series.
    inline const T *end() const;

    /// Get coefficient.
    //  Return a reference to coefficient of order [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    T &operator[](int n);

    /// Get coefficient.
    //  Return a reference to constant coefficient of order [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    const T &operator[](int n) const;


    /// Change sign of series.
    Taylor operator-() const;

    /// Multiply by scalar and assign.
    Taylor &operator*=(const T &);

    /// Divide by scalar and assign.
    Taylor &operator/=(const T &);

    /// Add series and assign.
    Taylor &operator+=(const Taylor &);

    /// Subtract series and assign.
    Taylor &operator-=(const Taylor &);


    /// Clear all coefficients.
    void clear();

    /// Integrate with respect to the variable.
    Taylor integrate() const;

    /// Return sum of series.
    T sum() const;


    /// Return order of this series.
    inline int getOrder() const;

private:

    // The representation of the coefficients.
    Array1D<T> itsCoeffs;
};


// Global functions.
// ------------------------------------------------------------------------

/// Add.
template <class T>
Taylor<T>
operator+(const Taylor<T> &, const Taylor<T> &);

/// Subtract.
template <class T>
Taylor<T>
operator-(const Taylor<T> &, const Taylor<T> &);

/// Multiply by scalar.
template <class T>
Taylor<T>
operator*(const Taylor<T> &, double);

/// Multiply by scalar.
template <class T>
Taylor<T>
operator*(double, const Taylor<T> &);

/// Divide by scalar.
template <class T>
Taylor<T>
operator/(const Taylor<T> &, double);

/// Poisson bracket of two Taylor seriess.
//  For this function the coefficients must have a Poisson bracket operation.
template <class T>
Taylor<T>
PoissonBracket(const Taylor<T> &, const Taylor<T> &);

/// Output.
template <class T>
std::ostream &operator<<(std::ostream &, const Taylor<T> &);

#include "FixedAlgebra/Taylor.hpp"

#endif // CLASSIC_Taylor_HH
