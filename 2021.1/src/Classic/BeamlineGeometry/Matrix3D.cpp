// ------------------------------------------------------------------------
// $RCSfile: Matrix3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Matrix3D
//   Matrix3D represents an 3-dimension matrix.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Matrix3D.h"
#include "BeamlineGeometry/Vector3D.h"


// Class Matrix3D.
// ------------------------------------------------------------------------

Matrix3D::Matrix3D(const Vector3D &a, const Vector3D &b, const Vector3D &c) {
    m[0][0] = a(0);
    m[0][1] = b(0);
    m[0][2] = c(0);
    m[1][0] = a(1);
    m[1][1] = b(1);
    m[1][2] = c(1);
    m[2][0] = a(2);
    m[2][1] = b(2);
    m[2][2] = c(2);
}


Matrix3D::Matrix3D(double x11, double x12, double x13,
                   double x21, double x22, double x23,
                   double x31, double x32, double x33) {
    m[0][0] = x11;
    m[0][1] = x12;
    m[0][2] = x13;
    m[1][0] = x21;
    m[1][1] = x22;
    m[1][2] = x23;
    m[2][0] = x31;
    m[2][1] = x32;
    m[2][2] = x33;
}


bool Matrix3D::operator==(const Matrix3D &rhs) const {
    return (m[0][0] == rhs(0, 0)) && (m[0][1] == rhs(0, 1)) && (m[0][2] == rhs(0, 2))
           && (m[1][0] == rhs(1, 0)) && (m[1][1] == rhs(1, 1)) && (m[1][2] == rhs(1, 2))
           && (m[2][0] == rhs(2, 0)) && (m[2][1] == rhs(2, 1)) && (m[2][2] == rhs(2, 2));
}


bool Matrix3D::operator!=(const Matrix3D &rhs) const {
    return (m[0][0] != rhs(0, 0)) || (m[0][1] != rhs(0, 1)) || (m[0][2] != rhs(0, 2))
           || (m[1][0] != rhs(1, 0)) || (m[1][1] != rhs(1, 1)) || (m[1][2] != rhs(1, 2))
           || (m[2][0] != rhs(2, 0)) || (m[2][1] != rhs(2, 1)) || (m[2][2] != rhs(2, 2));
}


Matrix3D &Matrix3D::operator+=(const Matrix3D &rhs) {
    m[0][0] += rhs.m[0][0];
    m[0][1] += rhs.m[0][1];
    m[0][2] += rhs.m[0][2];
    m[1][0] += rhs.m[1][0];
    m[1][1] += rhs.m[1][1];
    m[1][2] += rhs.m[1][2];
    m[2][0] += rhs.m[2][0];
    m[2][1] += rhs.m[2][1];
    m[2][2] += rhs.m[2][2];
    return *this;
}


Matrix3D &Matrix3D::operator-=(const Matrix3D &rhs) {
    m[0][0] -= rhs.m[0][0];
    m[0][1] -= rhs.m[0][1];
    m[0][2] -= rhs.m[0][2];
    m[1][0] -= rhs.m[1][0];
    m[1][1] -= rhs.m[1][1];
    m[1][2] -= rhs.m[1][2];
    m[2][0] -= rhs.m[2][0];
    m[2][1] -= rhs.m[2][1];
    m[2][2] -= rhs.m[2][2];
    return *this;
}


Matrix3D &Matrix3D::operator*=(const Matrix3D &rhs) {
    return (*this = *this * rhs);
}


Matrix3D &Matrix3D::operator*=(double rhs) {
    m[0][0] *= rhs;
    m[0][1] *= rhs;
    m[0][2] *= rhs;
    m[1][0] *= rhs;
    m[1][1] *= rhs;
    m[1][2] *= rhs;
    m[2][0] *= rhs;
    m[2][1] *= rhs;
    m[2][2] *= rhs;
    return *this;
}


Matrix3D Matrix3D::Identity() {
    return Matrix3D();
}


bool Matrix3D::isIdentity() const {
    return (m[0][0] == 1.0) && (m[0][1] == 0.0) && (m[0][2] == 0.0)
           && (m[1][0] == 0.0) && (m[1][1] == 1.0) && (m[1][2] == 0.0)
           && (m[2][0] == 0.0) && (m[2][1] == 0.0) && (m[2][2] == 1.0);
}


Matrix3D Matrix3D::inverse() const {
    double det =
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) +
        m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2]) +
        m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    return Matrix3D((m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
                    (m[2][1] * m[0][2] - m[2][2] * m[0][1]) / det,
                    (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,

                    (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
                    (m[2][2] * m[0][0] - m[2][0] * m[0][2]) / det,
                    (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,

                    (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
                    (m[2][0] * m[0][1] - m[2][1] * m[0][0]) / det,
                    (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det);
}


Matrix3D Matrix3D::transpose() const {
    Matrix3D t;

    for(int i = 0; i < 3; i++) {
        t(i, 0) = m[0][i];
        t(i, 1) = m[1][i];
        t(i, 2) = m[2][i];
    }

    return t;
}


Matrix3D operator+(const Matrix3D &lhs, const Matrix3D &rhs) {
    Matrix3D result(lhs);
    return result += rhs;
}


Matrix3D operator-(const Matrix3D &lhs, const Matrix3D &rhs) {
    Matrix3D result(lhs);
    return result += rhs;
}


Matrix3D operator*(const Matrix3D &lhs, const Matrix3D &rhs) {
    return Matrix3D(lhs(0, 0) * rhs(0, 0) + lhs(0, 1) * rhs(1, 0) + lhs(0, 2) * rhs(2, 0),
                    lhs(0, 0) * rhs(0, 1) + lhs(0, 1) * rhs(1, 1) + lhs(0, 2) * rhs(2, 1),
                    lhs(0, 0) * rhs(0, 2) + lhs(0, 1) * rhs(1, 2) + lhs(0, 2) * rhs(2, 2),

                    lhs(1, 0) * rhs(0, 0) + lhs(1, 1) * rhs(1, 0) + lhs(1, 2) * rhs(2, 0),
                    lhs(1, 0) * rhs(0, 1) + lhs(1, 1) * rhs(1, 1) + lhs(1, 2) * rhs(2, 1),
                    lhs(1, 0) * rhs(0, 2) + lhs(1, 1) * rhs(1, 2) + lhs(1, 2) * rhs(2, 2),

                    lhs(2, 0) * rhs(0, 0) + lhs(2, 1) * rhs(1, 0) + lhs(2, 2) * rhs(2, 0),
                    lhs(2, 0) * rhs(0, 1) + lhs(2, 1) * rhs(1, 1) + lhs(2, 2) * rhs(2, 1),
                    lhs(2, 0) * rhs(0, 2) + lhs(2, 1) * rhs(1, 2) + lhs(2, 2) * rhs(2, 2));
}


Vector3D operator*(const Matrix3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs(0, 0) * rhs(0) + lhs(0, 1) * rhs(1) + lhs(0, 2) * rhs(2),
                    lhs(1, 0) * rhs(0) + lhs(1, 1) * rhs(1) + lhs(1, 2) * rhs(2),
                    lhs(2, 0) * rhs(0) + lhs(2, 1) * rhs(1) + lhs(2, 2) * rhs(2));
}
