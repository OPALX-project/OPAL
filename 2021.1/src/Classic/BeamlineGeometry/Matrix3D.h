#ifndef CLASSIC_Matrix3D_HH
#define CLASSIC_Matrix3D_HH

// ------------------------------------------------------------------------
// $RCSfile: Matrix3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Matrix3D
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

class Vector3D;


// Class Matrix3D
// ------------------------------------------------------------------------
/// 3-dimensional matrix.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons they
//  are not implemented.

class Matrix3D {

public:

    /// Default constructor.
    //  Constructs identity.
    Matrix3D();

    /// Constructor.
    //  Use the three vectors (a,b,c) as column vectors.
    Matrix3D(const Vector3D &a, const Vector3D &b, const Vector3D &c);

    /// Constructor.
    //  Use the elements as matrix elements by rows.
    Matrix3D(double x11, double x12, double x13,
             double x21, double x22, double x23,
             double x31, double x32, double x33);

    bool operator==(const Matrix3D &) const;
    bool operator!=(const Matrix3D &) const;

    /// Get element.
    //  Return reference to matrix element (i,k).
    double &operator()(int i, int k)      ;

    /// Get element.
    //  Return value of matrix element (i,k).
    double operator()(int i, int k) const;

    /// Add and assign.
    Matrix3D &operator+=(const Matrix3D &rhs);

    /// Subtract and assign.
    Matrix3D &operator-=(const Matrix3D &rhs);

    /// Multiply and assign.
    Matrix3D &operator*=(const Matrix3D &rhs);

    /// Multiply and assign.
    Matrix3D &operator*=(double factor)      ;

    /// Make identity.
    static Matrix3D Identity();

    /// Inverse.
    Matrix3D inverse() const;

    /// Test for identity.
    //  Return true, if [b]this[/b] is an identity matrix.
    bool isIdentity() const;

    /// Transpose.
    Matrix3D transpose() const;

private:

    // matrix elements
    double m[3][3];
};


// Global operators.
// ------------------------------------------------------------------------

/// Add.
Matrix3D operator+(const Matrix3D &lhs, const Matrix3D &rhs);

/// Subtract.
Matrix3D operator-(const Matrix3D &lhs, const Matrix3D &rhs);

/// Multiply.
Matrix3D operator*(const Matrix3D &lhs, const Matrix3D &rhs);

/// Multiply.
Vector3D operator*(const Matrix3D &lhs, const Vector3D &rhs);


// Inline functions.
// ------------------------------------------------------------------------

inline Matrix3D::Matrix3D() {
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
}


inline double &Matrix3D::operator()(int i, int k)
{ return m[i][k]; }


inline double Matrix3D::operator()(int i, int k) const
{ return m[i][k]; }

#endif // __Matrix3D_HH
