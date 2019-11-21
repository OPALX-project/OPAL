#ifndef CLASSIC_FVector_HH
#define CLASSIC_FVector_HH

// ------------------------------------------------------------------------
// $RCSfile: FVector.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FVector<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FArray1D.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <functional>

// Template class FVector<T,N>
// ------------------------------------------------------------------------
/// A templated representation for vectors.
//  This class implements the arithmetic operations.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons
//  they are not implemented.

template<class T, int N>
class FVector: public FArray1D<T, N> {

public:

    /// Constructor.
    //  Construct zero FVector.
    FVector();

    /// Constructor.
    //  Set all vector elements to [b]t[/b].
    explicit FVector(const T &t);

    /// Constructor.
    //  Fill all vector elements from the C-array [b]t[/b].
    explicit FVector(const T *t);

    /// Conversion from one-dimensional array.
    FVector(const FArray1D<T, N> &);

    /// Change sign.
    FVector operator-() const;

    /// Multiply by scalar and assign.
    FVector &operator*=(const T &);

    /// Divide by scalar and assign.
    FVector &operator/=(const T &);

    /// Add FVector and assign.
    FVector &operator+=(const FVector &);

    /// Subtract FVector and assign.
    FVector &operator-=(const FVector &);
};


// ------------------------------------------------------------------------

/// Add.
template<class T, int N>
FVector<T, N> operator+(const FVector<T, N> &, const FVector<T, N> &);

/// Subtract.
template<class T, int N>
FVector<T, N> operator-(const FVector<T, N> &, const FVector<T, N> &);

/// Dot product.
template<class T, int N>
T operator*(const FVector<T, N> &, const FVector<T, N> &);

/// Multiply.
template<class T, int N>
FVector<T, N> operator*(const FVector<T, N> &, const T &);

/// Divide.
template<class T, int N>
FVector<T, N> operator/(const FVector<T, N> &, const T &);

/// Multiply.
template<class T, int N>
FVector<T, N> operator*(const T &, const FVector<T, N> &);

/// Euclidean norm.
template<class T, int N>
T euclidean_norm(const FVector<T, N> &);

/// Euclidean norm of diagonal matrix D times FVector V.
template<class T, int N>
T scaled_norm(const FArray1D<T, N> D, const FVector<T, N> &V);


// Implementation of template class FVector<T,N>.
// ------------------------------------------------------------------------

template<class T, int N>
FVector<T, N>::FVector():
    FArray1D<T, N>()
{}


template<class T, int N>
FVector<T, N>::FVector(const T &val):
    FArray1D<T, N>(val)
{}


template<class T, int N>
FVector<T, N>::FVector(const T *rhs):
    FArray1D<T, N>(rhs)
{}


template<class T, int N>
FVector<T, N>::FVector(const FArray1D<T, N> &rhs):
    FArray1D<T, N>(rhs)
{}


template<class T, int N>
FVector<T, N> FVector<T, N>::operator-() const {
    FVector<T, N> result;
    std::transform(this->begin(), this->end(), result.begin(), std::negate<T>());
    return result;
}


template<class T, int N>
FVector<T, N> &FVector<T, N>::operator*=(const T &val) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, val));
    return *this;
}


template<class T, int N>
FVector<T, N> &FVector<T, N>::operator/=(const T &val) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind(std::divides<T>(), std::placeholders::_1, val));
    return *this;
}


template<class T, int N>
FVector<T, N> &FVector<T, N>::operator+=(const FVector<T, N> &rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(),
                   std::plus<T>());
    return *this;
}


template<class T, int N>
FVector<T, N> &FVector<T, N>::operator-=(const FVector<T, N> &rhs) {
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(),
                   std::minus<T>());
    return *this;
}


// Global template operators.
// ------------------------------------------------------------------------

template<class T, int N>
FVector<T, N> operator+(const FVector<T, N> &lhs, const FVector<T, N> &rhs) {
    FVector<T, N> result;
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   result.begin(), std::plus<T>());
    return result;
}


template<class T, int N>
FVector<T, N> operator-(const FVector<T, N> &lhs, const FVector<T, N> &rhs) {
    FVector<T, N> result;
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   result.begin(), std::minus<T>());
    return result;
}


// NOTE: this is the standard FVector dot product,
// NOT memberwise multiplication.
template<class T, int N>
T operator*(const FVector<T, N> &lhs, const FVector<T, N> &rhs) {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), T(0));
}


template<class T, int N>
FVector<T, N> operator*(const FVector<T, N> &lhs, const T &x) {
    FVector<T, N> result(lhs);
    return result *= x;
}


template<class T, int N>
FVector<T, N> operator/(const FVector<T, N> &lhs, const T &x) {
    FVector<T, N> result(lhs);
    return result /= x;
}


// NOTE: this function assumes that multiplication for T is commutative.
template<class T, int N>
FVector<T, N> operator*(const T &x, const FVector<T, N> &lhs) {
    FVector<T, N> result(lhs);
    std::transform(lhs.begin(), lhs.end(), result.begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, x));
    return result;
}


template<class T, int N>
T euclidean_norm(const FVector<T, N> &V) {
    return sqrt(std::inner_product(V.begin(), V.end(), V.begin(), T(0)));
}


template<class T, int N>
T scaled_norm(const FArray1D<T, N> D, const FVector<T, N> &V) {
    T sum(0);

    for(int i = 0; i < V.size(); ++i) {
        double dv = D[i] * V[i];
        sum += dv * dv;
    }

    return sqrt(sum);
}

#endif //  CLASSIC_FVector_HH
