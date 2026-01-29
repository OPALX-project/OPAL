#ifndef CLASSIC_EMField_HH
#define CLASSIC_EMField_HH

// ------------------------------------------------------------------------
// $RCSfile: EMField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Classes:
//   Point3D:      a point in 3-dimensional space.
//   EVector:      an electric field vector.
//   BVector:      a magnetic field vector.
//   EBVectors:    an electromagnetic field, containing E and B.
//   EMField:      a virtual base class for electromagnetic fields.
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:36 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class Point3D
// ------------------------------------------------------------------------
/// A point in 3 dimensions.

class Point3D {

public:

    /// Constructor.
    //  Construct the point from its coordinates (x,y,z) in m.
    Point3D(double x, double y, double z);

    /// Return coordinate.
    //  Return the x-coordinate for the point in m.
    double getX() const;

    /// Return coordinate.
    //  Return the y-coordinate for the point in m.
    double getY() const;

    /// Return coordinate.
    //  Return the z-coordinate for the point in m.
    double getZ() const;

private:

    double x, y, z;
};


// Class EVector
// ------------------------------------------------------------------------
/// An electric field vector.

class EVector {

public:

    /// Constructor.
    //  Construct the field vector from is components (Ex,Ey,Ez) in A/m.
    EVector(double, double, double);

    /// Scale.
    //  Multiply the field vector by [b]scalar[/b].
    EVector operator*(double scalar) const;

    /// Get component.
    //  Return the x-component of the field in A/m.
    double getEx() const;

    /// Get component.
    //  Return the y-component of the field in A/m.
    double getEy() const;

    /// Get component.
    //  Return the z-component of the field in A/m.
    double getEz() const;

private:

    double Ex, Ey, Ez;
};


// Class BVector
// ------------------------------------------------------------------------
/// A magnetic field vector.

class BVector {

public:

    /// Constructor.
    //  Construct the field vector from its components (Bx,By,Bz) in T.
    BVector(double, double, double);

    /// Scale.
    //  Multiply the field vector by [b]scalar[/b].
    BVector operator*(double scalar) const;

    /// Get component.
    //  Return the x-component of the field in T.
    double getBx() const;

    /// Get component.
    //  Return the y-component of the field in T.
    double getBy() const;

    /// Get component.
    //  Return the z-component of the field in T.
    double getBz() const;

private:

    double Bx, By, Bz;
};


// Class EBVectors
// ------------------------------------------------------------------------
/// A representation of an electromagnetic field.
//  Contains both an electric and a magnetic field vector.
//  These represent the instantaneous electromagnetic field in a given point.

class EBVectors {

public:

    /// Constructor.
    //  Construct the field from the two field vectors [b]E[/b] (in A/m)
    //  and [b]B[/b] (in T).
    EBVectors(const EVector &E, const BVector &B);

    /// Get component.
    //  Return the electric field vector in A/m.
    EVector getE() const;

    /// Get component.
    //  Return the x-component of the electric field in A/m.
    double getEx() const;

    /// Get component.
    //  Return the y-component of the electric field in A/m.
    double getEy() const;

    /// Get component.
    //  Return the z-component of the electric field in A/m.
    double getEz() const;

    /// Get field.
    //  Return the magnetic field vector in T.
    BVector getB() const;

    /// Get component.
    //  Return the x-component of the magnetic field in T.
    double getBx() const;

    /// Get component.
    //  Return the y-component of the magnetic field in T.
    double getBy() const;

    /// Get component.
    //  Return the z-component of the magnetic field in T.
    double getBz() const;

private:

    EVector E;
    BVector B;
};


// Class EMField
// ------------------------------------------------------------------------
/// Abstract base class for electromagnetic fields.
//  This class represent a time-varying position dependent electromagnetic
//  field.  Derived classes include time-constant fields as well as
//  spatially homogeneous fields.

class EMField {

public:

    /// Default constructor.
    //  Construct zero field.
    EMField();

    EMField(const EMField &right);
    virtual ~EMField();
    const EMField &operator=(const EMField &right);

    /// Get field.
    //  Return the time-independent part of the electric field in point [b]P[/b].
    //  This default action returns a zero field.
    virtual EVector Efield(const Point3D &P) const;

    /// Get field.
    //  Return the time-independent part of the magnetic field in point [b]P[/b].
    //  This default action returns a zero field.
    virtual BVector Bfield(const Point3D &P) const;

    /// Get field.
    //  Return the electric field at time [b]t[/b] in point [b]P[/b].
    //  This default action returns the static part  EField(P).
    virtual EVector Efield(const Point3D &P, double t) const;

    /// Get field.
    //  Return the magnetic field at time [b]t[/b] in point [b]P[/b].
    //  This default action returns the static part  BField(P).
    virtual BVector Bfield(const Point3D &P, double t) const;

    /// Get field.
    //  Return the static part of the field pair (E,B) in point [b]P[/b].
    //  This default action returns  EBVectors(EField(P), BField(P)).
    virtual EBVectors EBfield(const Point3D &P) const;

    /// Get field.
    //  Return the field pair (E,B) at time [b]t[/b] in point [b]P[/b].
    //  This default action returns the static part  EBfield(P).
    virtual EBVectors EBfield(const Point3D &P, double t) const;

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    //  This method must be defined for each derived class.
    virtual void scale(double scalar) = 0;

    /// The constant representing a zero electric field.
    static const EVector ZeroEfield;

    /// The constant representing a zero magnetic field.
    static const BVector ZeroBfield;

    /// The constant representing a zero electromagnetic field.
    static const EBVectors ZeroEBfield;
};

#endif // CLASSIC_EMField_HH
