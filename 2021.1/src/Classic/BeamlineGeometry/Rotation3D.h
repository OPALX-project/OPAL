#ifndef CLASSIC_Rotation3D_HH
#define CLASSIC_Rotation3D_HH

// ------------------------------------------------------------------------
// $RCSfile: Rotation3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Rotation3D
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


// Class Rotation3D
// -----------------------------------------------------------------------
/// Rotation in 3-dimensional space.
//  There are two possible external representations:
//  [DL]
//  [DT]Matrix3D R[DD]
//  The matrix R is an element of SO(3). Its column vectors define the
//  unit vectors of the rotated coordinate system, expressed in the
//  original system. This representation is used internally and can only
//  be read from a Rotation3D object. To build such an object, use the
//  other representation.
//  [DT]Vector3D V[DD]
//  The direction of V defines the axis of rotation, and its length the
//  sine of half the rotation angle. A zero vector represents the identity.
//  [/DL]
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons they
//  are not implemented.

class Rotation3D {

public:

    /// Default constructor.
    //  Constructs identity.
    Rotation3D();

    /// Constructor.
    //  Use axis vector (vx,vy,vz).
    Rotation3D(double vx, double vy, double vz);

    bool operator==(const Rotation3D &) const;
    bool operator!=(const Rotation3D &) const;

    /// Multiply and assign.
    Rotation3D &operator*=(const Rotation3D &rhs);

    /// Multiply.
    Rotation3D operator*(const Rotation3D &rhs) const;

    /// Rotate.
    //  Return vector [b]rhs[/b] rotated by the rotation [b]this[/b].
    Vector3D operator*(const Vector3D &rhs) const;

    /// Get component.
    //  Return element (i,k) of the orthogonal matrix.
    double operator()(int i, int k) const;

    /// Get axis vector.
    //  Return the axis vector.
    Vector3D getAxis() const;

    /// Get axis vector.
    //  Return the axis vector in (vx,vy,vz).
    void getAxis(double &vx, double &vy, double &vz) const;

    /// Inversion.
    Rotation3D inverse() const;

    /// Test for identity.
    bool isIdentity() const;

    /// Test for rotation.
    //  Return true, if [b]this[/b] is a pure x-rotation.
    bool isPureXRotation() const;

    /// Test for rotation.
    //  Return true, if [b]this[/b] is a pure y-rotation.
    bool isPureYRotation() const;

    /// Test for rotation.
    //  Return true, if [b]this[/b] is a pure z-rotation.
    bool isPureZRotation() const;

    /// Make identity.
    static Rotation3D Identity();

    /// Make rotation.
    //  Construct pure rotation around x-axis.
    static Rotation3D XRotation(double angle);

    /// Make rotation.
    //  Construct pure rotation around y-axis.
    static Rotation3D YRotation(double angle);

    /// Make rotation.
    //  Construct pure rotation around z-axis.
    static Rotation3D ZRotation(double angle);

private:

    // Representation.
    Matrix3D R;
};


// Inline functions.
// ------------------------------------------------------------------------

inline Rotation3D::Rotation3D():
    R()
{}


inline bool Rotation3D::operator==(const Rotation3D &rhs) const
{ return (R == rhs.R); }


inline bool Rotation3D::operator!=(const Rotation3D &rhs) const
{ return (R != rhs.R); }


inline double Rotation3D::operator()(int row, int col) const
{ return R(row, col); }


inline Rotation3D &Rotation3D::operator*=(const Rotation3D &rhs) {
    R *= rhs.R;
    return *this;
}


inline Rotation3D Rotation3D::operator*(const Rotation3D &rhs) const {
    Rotation3D result(*this);
    return result *= rhs;
}


inline Vector3D Rotation3D::operator*(const Vector3D &rhs) const {
    return R * rhs;
}


inline bool Rotation3D::isIdentity() const {
    return R.isIdentity();
}

#endif // CLASSIC_Rotation3D_HH
