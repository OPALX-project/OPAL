#ifndef CLASSIC_Vector_HH
#define CLASSIC_Vector_HH

// ------------------------------------------------------------------------
// $RCSfile: Vector.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.4.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Vector
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:50:59 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <functional>

// Template class Vector<T>
// ------------------------------------------------------------------------
/// Vector.
//  A templated representation for vectors.
//  This class implements the arithmetic operations.

template<class T>
class Vector: public Array1D<T> {

public:

    /// Default constructor.
    //  Constructs vector of zero length.
    Vector();

    /// Constructor.
    //  Reserve [b]n[/b] members and leave them undefined.
    explicit Vector(int n);

    /// Constructor.
    //  Reserve [b]n[/b] members and set them to [b]t[/b].
    Vector(int n, const T &t);

    /// Conversion.
    Vector(const Array1D<T> &);

    Vector(const Vector<T> &);
    ~Vector();
    Vector<T> &operator=(const Vector<T> &);

    /// Change sign of vector.
    Vector<T> operator-() const;

    /// Multiply by scalar and assign.
    Vector<T> &operator*=(const T &);

    /// Divide by scalar and assign.
    Vector<T> &operator/=(const T &);

    /// Add vector and assign.
    Vector<T> &operator+=(const Vector<T> &);

    /// Subtract vector and assign.
    Vector<T> &operator-=(const Vector<T> &);
};


// Global Operators and Functions on Vector<T>
// ------------------------------------------------------------------------

/// Vector addition.
template<class T> Vector<T> operator+(const Vector<T> &, const Vector<T> &);

/// Vector subtraction.
template<class T> Vector<T> operator-(const Vector<T> &, const Vector<T> &);

/// Vector dot product.
template<class T> T operator*(const Vector<T> &, const Vector<T> &);

/// Vector multiplied by scalar.
template<class T> Vector<T> operator*(const Vector<T> &, const T &);

/// Vector divided by scalar.
template<class T> Vector<T> operator/(const Vector<T> &, const T &);

/// Scalar multiplied by vector.
template<class T> Vector<T> operator*(const T &, const Vector<T> &);

/// Euclidean norm.
template<class T> T euclidean_norm(const Vector<T> &);

/// Euclidean norm of diagonal matrix D times vector V.
template<class T> T scaled_norm(const Array1D<T> D, const Vector<T> &V);


// Implementation of template class Vector<T>.
// ------------------------------------------------------------------------

template<class T>
Vector<T>::Vector():
    Array1D<T>()
{}


template<class T>
Vector<T>::Vector(int n):
    Array1D<T>(n)
{}


template<class T>
Vector<T>::Vector(int n, const T &val):
    Array1D<T>(n, val)
{}


template<class T>
Vector<T>::Vector(const Vector<T> &array):
    Array1D<T>(array)
{}


template<class T>
Vector<T>::Vector(const Array1D<T> &array):
    Array1D<T>(array)
{}


template<class T>
Vector<T>::~Vector()
{}


template<class T>
Vector<T> &Vector<T>::operator=(const Vector<T> &right) {
    Array1D<T>::operator=(right);
    return *this;
}


template<class T>
Vector<T> Vector<T>::operator-() const {
    Vector<T> result(this->size());
    std::transform(this->begin(), this->end(), result.begin(), std::negate<T>());
    return result;
}


template<class T>
Vector<T> &Vector<T>::operator*=(const T &val) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, val));
    return *this;
}


template<class T>
Vector<T> &Vector<T>::operator/=(const T &val) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind(std::divides<T>(), std::placeholders::_1, val));
    return *this;
}


template<class T>
Vector<T> &Vector<T>::operator+=(const Vector<T> &array) {
    std::transform(this->begin(), this->end(), array.begin(), this->begin(),
                   std::plus<T>());
    return *this;
}


template<class T>
Vector<T> &Vector<T>::operator-=(const Vector<T> &array) {
    std::transform(this->begin(), this->end(), array.begin(), this->begin(),
                   std::minus<T>());
    return *this;
}


// Global template operators.
// ------------------------------------------------------------------------

template<class T>
Vector<T> operator+(const Vector<T> &V1, const Vector<T> &V2) {
    Vector<T> result(V1.size());
    std::transform(V1.begin(), V1.end(), V2.begin(), result.begin(),
                   std::plus<T>());
    return result;
}


template<class T>
Vector<T> operator-(const Vector<T> &V1, const Vector<T> &V2) {
    Vector<T> result(V1.size());
    std::transform(V1.begin(), V1.end(), V2.begin(), result.begin(),
                   std::minus<T>());
    return result;
}


// NOTE: this is the standard vector dot product,
// NOT memberwise multiplication.
template<class T>
T operator*(const Vector<T> &V1, const Vector<T> &V2) {
    return std::inner_product(V1.begin(), V1.end(), V2.begin(), T(0));
}


template<class T>
Vector<T> operator*(const Vector<T> &V1, const T &x) {
    Vector<T> result = V1;
    return result *= x;
}


template<class T>
Vector<T> operator/(const Vector<T> &V1, const T &x) {
    Vector<T> result = V1;
    return result /= x;
}


// NOTE: this function assumes that multiplication for T is commutative.
template<class T>
Vector<T> operator*(const T &x, const Vector<T> &V1) {
    Vector<T> result = V1;
    std::transform(V1.begin(), V1.end(), result.begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, x));
    return result;
}


template<class T> T euclidean_norm(const Vector<T> &V) {
    return sqrt(std::inner_product(V.begin(), V.end(), V.begin(), T(0)));
}


template<class T> T scaled_norm(const Array1D<T> D, const Vector<T> &V) {
    if(D.size() != V.size()) {
        throw SizeError("scaled_norm()", "Inconsistent dimensions.");
    }

    T sum(0);
    for(int i = 0; i < V.size(); ++i) {
        double dv = D[i] * V[i];
        sum += dv * dv;
    }
    return std::sqrt(sum);
}

#endif //  CLASSIC_Vector_HH
