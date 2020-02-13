#ifndef CLASSIC_Matrix_HH
#define CLASSIC_Matrix_HH

// ------------------------------------------------------------------------
// $RCSfile: Matrix.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Matrix
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:05 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algebra/Array2D.h"
#include "Algebra/Vector.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <numeric>


// Template class Matrix<T>
// ------------------------------------------------------------------------
/// Matrix.
//  A templated representation for matrices.
//  This class implements the basic algebraic operations on general matrices,
//  which need not be square.

template<class T>
class Matrix: public Array2D<T> {

public:

    /// Default constructor.
    //  Constructs undefined matrix.
    Matrix();

    /// Conversion.
    //  From two-dimensional array.
    Matrix(const Array2D<T> &);

    /// Constructor.
    //  Reserve [b]rows[/b] x [b]cols[/b] elements and leave them undefined.
    Matrix(int rows, int cols);

    /// Constructor.
    //  Reserve [b]rows[/b] x [b]cols[/b] elements and set them to t.
    Matrix(int rows, int cols, const T &t);

    Matrix(const Matrix &);
    ~Matrix();
    Matrix<T> &operator=(const Matrix<T> &);

    /// Convert and assign.
    Matrix<T> &operator=(const Array2D<T> &);

    /// Multiply by scalar and assign.
    Matrix<T> &operator*=(const T &);

    /// Divide by scalar and assign.
    Matrix<T> &operator/=(const T &);

    /// Add matrix and assign..
    Matrix<T> &operator+=(const Matrix<T> &);

    /// Subtract matrix and assign.
    Matrix<T> &operator-=(const Matrix<T> &);

    /// Add scalar.
    //  Add the unit matrix times [b]t[/b] and assign.
    //  Throw SizeError, if matrix is not square.
    Matrix<T> &operator+=(const T &);

    /// Subtract scalar.
    //  Subtract the unit matrix times [b]t[/b] and assign.
    //  Throw SizeError, if matrix is not square.
    Matrix<T> &operator-=(const T &);

    /// Multiply by matrix and assign.
    Matrix<T> &operator*=(const Matrix<T> &);

    /// Matrix transpose.
    Matrix<T> transpose() const;

    /// Matrix product.
    Matrix<T> dotm(const Matrix<T> &) const;

    /// Matrix times column vector.
    Vector<T> dotcv(const Vector<T> &) const;

    /// Row vector times matrix.
    Vector<T> dotrv(const Vector<T> &) const;
};


// Global Operators on Matrix<T>
// ------------------------------------------------------------------------

/// Matrix addition.
template <class T> Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);

/// Matrix subtraction.
template <class T> Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);

/// Matrix unary minus.
template <class T> Matrix<T> operator-(const Matrix<T> &);

/// Matrix multiply.
template <class T> Matrix<T> operator*(const Matrix<T> &, const Matrix<T> &);

/// Matrix times column vector.
template <class T> Vector<T> operator*(const Matrix<T> &, const Vector<T> &);

/// Row vector times matrix.
template <class T> Vector<T> operator*(const Vector<T> &, const Matrix<T> &);

/// Matrix times scalar.
template <class T> Matrix<T> operator*(const Matrix<T> &, const T &);

/// Scalar times matrix.
template <class T> Matrix<T> operator*(const T &, const Matrix<T> &);

/// Matrix divided by scalar.
template <class T> Matrix<T> operator/(const Matrix<T> &, const T &);


// Public methods of Matrix<T>.
// ------------------------------------------------------------------------

template<class T>
Matrix<T>::Matrix():
    Array2D<T>()
{}


template<class T>
Matrix<T>::Matrix(const Matrix<T> &rhs):
    Array2D<T>(rhs)
{}


template<class T>
Matrix<T>::Matrix(const Array2D<T> &rhs):
    Array2D<T>(rhs)
{}


template<class T>
Matrix<T>::Matrix(int rows, int cols):
    Array2D<T>(rows, cols)
{}


template<class T>
Matrix<T>::Matrix(int rows, int cols, const T &/*val*/):
    Array2D<T>(rows, cols, T(0)) {
    for(int i = std::min(rows, cols); i-- > 0;) {
        this->data[i * (cols + 1)] = T(1);
    }
}


template<class T>
Matrix<T>::~Matrix()
{}


template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs) {
    Array2D<T>::operator=(rhs);
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator=(const Array2D<T> &rhs) {
    Array2D<T>::operator=(rhs);
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator*=(const T &val) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind2nd(std::multiplies<T>(), val));
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator/=(const T &val) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind2nd(std::divides<T>(), val));
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &array) {
    if(this->nrows() != array.nrows()  ||  this->ncols() != array.ncols()) {
        throw SizeError("Matrix::operator+=()", "Dimensions inconsistent.");
    }
    std::transform(this->begin(), this->end(), array.begin(), this->begin(),
                   std::plus<T>());
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &array) {
    if(this->nrows() != array.nrows()  ||  this->ncols() != array.ncols()) {
        throw SizeError("Matrix::operator-=()", "Dimensions inconsistent.");
    }
    std::transform(this->begin(), this->end(), array.begin(), this->begin(),
                   std::minus<T>());
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator+=(const T &val) {
    if(this->nrows() != this->ncols()) {
        throw SizeError("Matrix::operator+=()", "Matrix is not square.");
    }
    int n = std::min(this->nrows(), this->ncols());
    for(int i = 0; i < n; i++) this->col_begin(i)[i] += val;
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator-=(const T &val) {
    if(this->nrows() != this->ncols()) {
        throw SizeError("Matrix::operator-=()", "Matrix is not square.");
    }
    int n = std::min(this->nrows(), this->ncols());
    for(int i = 0; i < n; i++) this->col_begin(i)[i] -= val;
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &array) {
    *this = (*this).dotm(array);
    return *this;
}


template <class T>
Matrix<T> Matrix<T>::transpose() const {
    int nr = this->nrows();
    int nc = this->ncols();
    Matrix<T> R(nc, nr);

    for(int i = 0; i < nr; ++i) {
        std::copy(this->row_begin(i), this->row_end(i), R.col_begin(i));
    }

    return R;
}


// Global functions acting on Matrix<T>.
// ------------------------------------------------------------------------

template<class T>
Matrix<T> operator+(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    Matrix<T> result = lhs;
    return result += rhs;
}


template<class T>
Matrix<T> operator-(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    Matrix<T> result = lhs;
    return result -= rhs;
}


template<class T>
Matrix<T> operator-(const Matrix<T> arg) {
    Matrix<T> result(arg.nrows(), arg.ncols());
    std::transform(arg.begin(), arg.end(), result.begin(), std::negate<T>());
    return result;
}


template<class T>
Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    return lhs.dotm(rhs);
}


template<class T>
Vector<T> operator*(const Matrix<T> &M, const Vector<T> &V) {
    return M.dotcv(V);
}


template<class T>
Vector<T> operator*(const Vector<T> &V, const Matrix<T> &M) {
    return M.dotrv(V);
}


template<class T>
Matrix<T> operator*(const Matrix<T> &M, const T &val) {
    Matrix<T> R(M);
    return R *= val;
}


template<class T>
Matrix<T> operator*(const T &val, const Matrix<T> &M) {
    Matrix<T> R(M);
    return R *= val;
}


template<class T>
Matrix<T> operator/(const Matrix<T> &M, const T &val) {
    Matrix<T> R(M);
    return R /= val;
}


// Protected methods of Matrix<T>.
// ------------------------------------------------------------------------

template<class T>
Matrix<T> Matrix<T>::dotm(const Matrix<T> &B) const {
    if(this->ncols() != B.nrows()) {
        throw SizeError("Matrix::dotm()", "Inconsistent dimensions.");
    }

    const int nr1 = this->nrows();
    const int nc2 = B.ncols();
    Matrix<T> R(nr1, nc2, T(0.0));

    for(int i = 0; i < nr1; i++) {
        typename Array2D<T>::const_row_iterator i1 = this->row_begin(i);
        typename Array2D<T>::const_row_iterator i2 = this->row_end(i);

        for(int j = 0; j < nc2; j++) {
            R(i, j) = std::inner_product(i1, i2, B.col_begin(j), T(0));
        }
    }

    return R;
}

template<class T>
Vector<T> Matrix<T>::dotcv(const Vector<T> &B) const {
    if(this->ncols() != B.size()) {
        throw SizeError("Matrix::dotcv()", "Inconsistent dimensions.");
    }

    int nr = this->nrows();
    Vector<T> R(nr, T(0));

    for(int i = 0; i < nr; ++i) {
        R(i) = std::inner_product(this->row_begin(i), this->row_end(i), B.begin(),
                                  T(0));
    }

    return R;
}

template<class T>
Vector<T> Matrix<T>::dotrv(const Vector<T> &B) const {
    if(this->nrows() != B.size()) {
        throw SizeError("Matrix::dotrv()", "Inconsistent dimensions.");
    }

    int nc = this->ncols();
    Vector<T> R(nc, T(0));

    for(int j = 0; j < nc; j++) {
        R(j) = std::inner_product(B.begin(), B.end(), this->col_begin(j), T(0));
    }

    return R;
}

#endif // CLASSIC_Matrix_HH
