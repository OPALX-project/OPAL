#ifndef CLASSIC_Vector3D_HH
#define CLASSIC_Vector3D_HH

// ------------------------------------------------------------------------
// $RCSfile: Vector3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Vector3D
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class Vector3D
// ------------------------------------------------------------------------
/// A 3-dimension vector.
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons they
//  are not implemented.

class Vector3D {

public:

    /// Default constructor.
    //  Construct null vector.
    Vector3D();

    /// Constructor.
    //  Use components (x,y,z).
    Vector3D(double x, double y, double z);

    bool operator==(const Vector3D &) const;
    bool operator!=(const Vector3D &) const;

    /// Add and assign.
    Vector3D &operator+=(const Vector3D &vector);

    /// Subtract and assign.
    Vector3D &operator-=(const Vector3D &vector);

    /// Scale and assign.
    Vector3D &operator*=(double factor);

    /// Negative vector.
    Vector3D operator-() const;

    /// Get component.
    //  Return a reference to component [b]i[/b].
    double &operator()(int i);

    /// Get component.
    //  Return the value of component [b]i[/b].
    double  operator()(int i) const;

    /// Get components.
    //  Return the components (x,y,z).
    void getComponents(double &x, double &y, double &z) const;

    /// Get component.
    //  Return the component [b]x[/b].
    double getX() const;

    /// Get component.
    //  Return the component [b]y[/b].
    double getY() const;

    /// Get component.
    //  Return the component [b]z[/b].
    double getZ() const;

    /// Set to zero.
    void clear();

    /// Test for zero.
    bool isZero() const;

    /// Set component.
    //  Assign the component [b]x[/b].
    void setX(double);

    /// Set component.
    //  Assign the component [b]y[/b].
    void setY(double);

    /// Set component.
    //  Assign the component [b]z[/b].
    void setZ(double);

protected:

    // Vector components.
    double v[3];
};


// External functions.
// ------------------------------------------------------------------------

/// Add.
extern Vector3D operator+(const Vector3D &a, const Vector3D &b);

/// Subtract.
extern Vector3D operator-(const Vector3D &a, const Vector3D &b);

/// Multiply.
extern Vector3D operator*(const Vector3D &a, double factor);

/// Multiply.
extern Vector3D operator*(double factor, const Vector3D &a);

/// Vector cross product.
extern Vector3D cross(const Vector3D &a, const Vector3D &b);

/// Vector dot product.
extern double dot(const Vector3D &a, const Vector3D &b);


// Inline functions.
// ------------------------------------------------------------------------

inline Vector3D::Vector3D()
{ v[0] = v[1] = v[2] = 0.0; }


inline double &Vector3D::operator()(int i)
{ return v[i]; }


inline double Vector3D::operator()(int i) const
{ return v[i]; }


inline double Vector3D::getX() const
{ return v[0]; }


inline double Vector3D::getY() const
{ return v[1]; }


inline double Vector3D::getZ() const
{ return v[2]; }


inline void Vector3D::setX(double x)
{ v[0] = x; }


inline void Vector3D::setY(double y)
{ v[1] = y; }


inline void Vector3D::setZ(double z)
{ v[2] = z; }

#endif // CLASSIC_Vector3D_HH
