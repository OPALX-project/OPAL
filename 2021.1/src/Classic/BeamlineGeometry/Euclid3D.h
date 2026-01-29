#ifndef CLASSIC_Euclid3D_HH
#define CLASSIC_Euclid3D_HH

// ------------------------------------------------------------------------
// $RCSfile: Euclid3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Euclid3D
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Rotation3D.h"
#include "BeamlineGeometry/Vector3D.h"


// Class Euclid3D
// ------------------------------------------------------------------------
/// Displacement and rotation in space.
//  An arbitrary 3-dimension translation and rotation.
//  The translation is a 3-dimensional vector, the rotation is represented
//  as a 3-dimensional orthonormal matrix.
//  {P}
//  If we think of an Euler3D object as an operator  E  which transforms a
//  point in the super-frame  X  to the local frame  X':
//  {P}
//  X' = E(X),
//  {P}
//  then we can define a dot product of two successive operations in the
//  following way: Given two Euclid3D objects (operators)  E1  and  E2,
//  we can define
//  {P}
//  E3 = E2 . E1
//  {P}
//  using
//  {P}
//  X''= E3(X) = E2(E1(X)).
//  {P}
//  The definition of a multiplicative inverse naturally follows as
//  {P}
//  E . inv(E) = I.
//  {P}
//  If  R  represents the 3-d rotation and  X  the vector specified by
//  (x,y,z),  then a transformation of a point  P  in the super-frame by
//  the Euclid3D is given by
//  {P}
//  P -> R . P + X.
//  {P}
//  For the representation of rotations, see Rotation3D.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons they
//  are not implemented.
//  An alternate representation for the rotation is given by a vector
//  pointing in the direction of the axis of rotation, whose sense is
//  given by the right-hand rule, and whose length is equal to the angle
//  of rotation.

class Euclid3D {
public:

    /// Default constructor.
    //  Construct identity.
    Euclid3D();

    /// Constructor/
    //  Use displacement [b]V[/b] and rotation [b]R[/b].
    Euclid3D(const Vector3D &V, const Rotation3D &R);

    /// Constructor.
    //  Use displacement (x,y,z) and rotation (vx,vy,vz).
    Euclid3D(double x, double y, double z, double vx, double vy, double vz);

    bool operator==(const Euclid3D &) const;
    bool operator!=(const Euclid3D &) const;

    /// Unpack.
    //  Return the displacement in (x,y,z) and the rotation axis in (vx,vy,vz).
    void getAll(double &x, double &y, double &z,
                double &vx, double &vy, double &vz) const;

    /// Get displacement.
    //  Return x-component of displacement.
    double getX() const;

    /// Get displacement.
    //  Return x-component of displacement.
    double getY() const;

    /// Get displacement.
    //  Return x-component of displacement.
    double getZ() const;

    /// Get displacement.
    //  Return displacement as a vector.
    const Vector3D &getVector() const;

    /// Get rotation.
    //  Return rotation as an orthogonal matrix.
    const Rotation3D &getRotation() const;

    /// Get component.
    //  Return the component (row,col) of the rotation matrix.
    double M(int row, int col) const;

    /// Set displacement.
    //  Assign x-component of displacement.
    void setX(double x);

    /// Set displacement.
    //  Assign y-component of displacement.
    void setY(double y);

    /// Set displacement.
    //  Assign z-component of displacement.
    void setZ(double z);

    /// Set displacement.
    //  Assign displacement as a vector.
    void setDisplacement(const Vector3D &V);

    /// Set rotation.
    //  Assign rotation as an orthogonal matrix.
    void setRotation(const Rotation3D &R);

    /// Dot product.
    //  Return composition of [b]this[/b] with [b]rhs[/b].
    Euclid3D dot(const Euclid3D &rhs) const;

    /// Dot product with assign.
    //  Return composition of [b]this[/b] with [b]rhs[/b].
    const Euclid3D &dotBy(const Euclid3D &rhs);

    /// Dot product.
    //  Return composition of [b]this[/b] with [b]rhs[/b].
    Euclid3D operator*(const Euclid3D &rhs) const;

    /// Dot product with assign.
    const Euclid3D &operator*=(const Euclid3D &rhs);

    /// Inverse.
    Euclid3D inverse() const;

    /// Test for identity.
    bool isIdentity() const;

    /// Test for translation.
    //  Return true, if [b]this[/b] is a pure translation.
    bool isPureTranslation() const;

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
    static Euclid3D identity();

    /// Make translation.
    // Use displacement (x,y,z).
    static Euclid3D translation(double x, double y, double z);

    /// Make rotation.
    //  Construct pure rotation by [b]angle[/b] around the x-axis;
    static Euclid3D XRotation(double angle);

    /// Make rotation.
    //  Construct pure rotation by [b]angle[/b] around the y-axis;
    static Euclid3D YRotation(double angle);

    /// Make rotation.
    //  Construct pure rotation by [b]angle[/b] around the z-axis;
    static Euclid3D ZRotation(double angle);

private:

    // Implementation.
    Vector3D V;                // The displacement.
    Rotation3D R;              // The rotation.
    mutable bool is_identity;  // True, when [b]this[/b] is identity.
};


// Global Function.
// ------------------------------------------------------------------------

/// Euclidean inverse.
inline Euclid3D Inverse(const Euclid3D &t) {
    return t.inverse();
}


// Trivial functions are inlined.
// ------------------------------------------------------------------------

inline Euclid3D::Euclid3D():
    V(), R(), is_identity(true)
{}

inline double Euclid3D::getX() const {
    return V.getX();
}

inline double Euclid3D::getY() const {
    return V.getY();
}

inline double Euclid3D::getZ() const {
    return V.getZ();
}

inline double Euclid3D::M(int row, int col) const {
    return R(row, col);
}


inline bool Euclid3D::isIdentity() const {
    return is_identity;
}

#endif // CLASSIC_Euclid3D_HH
