#ifndef CLASSIC_BMultipoleField_HH
#define CLASSIC_BMultipoleField_HH

// ------------------------------------------------------------------------
// $RCSfile: BMultipoleField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BMultipoleField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/StaticMagneticField.h"

// Class BMultipoleField
// ------------------------------------------------------------------------
/// The magnetic field of a multipole.
//  The static components are defined by the equation
//  {center}
//  By + i * Bx = sum (Bn + i * An) (x + i * y) ** (n - 1)
//  {/center}
//  where the index  n  runs from  1  (dipole) to  N  (2N-pole).
//  The field component n is stored in the pairs[n-1],
//  The array  pairs  has the dimension N.

class BMultipoleField: public StaticMagneticField {

public:

    /// Default constructor.
    //  Constructs a null field.
    BMultipoleField();

    BMultipoleField(const BMultipoleField &);
    virtual ~BMultipoleField();
    BMultipoleField &operator=(const BMultipoleField &);

    /// Get field.
    //  Return the time-independent part of the magnetic field in point [b]P[/b].
    //  This override forces implementation in derived classes.
    virtual BVector Bfield(const Point3D &P) const;

    /// Get field.
    //  Return the magnetic field at time [b]t[/b] in point [b]P[/b].
    //  This override forces implementation in derived classes.
    virtual BVector Bfield(const Point3D &P, double t) const;

    /// Get component.
    //  Return the value of the n'th normal multipole component in  T/m**(n-1).
    //  If [b]n[/b] is larger than N, the return value is zero.
    double getNormalComponent(int n) const;

    /// Get component.
    //  Return the value of the n'th skew multipole component in T/m**(n-1).
    //  If [b]n[/b] is larger than N, the return value is zero.
    double getSkewComponent(int n) const;

    /// Set component.
    //  Assign the value of the n'th normal multipole component in T/m**(n-1).
    //  If [b]n[/b] is larger than N, the new field component is inserted.
    void setNormalComponent(int n, double Bn);

    /// Set component.
    //  Assign the value of the n'th skew multipole component in T/m**(n-1).
    //  If [b]n[/b] is larger than N, the new field component is inserted.
    void setSkewComponent(int n, double Bn);

    /// Get component.
    //  Return the value of the n'th normal multipole component in T/m**(n-1).
    //  Fast version: [b]n[/b] is not checked.
    //  the result is undefined if it is larger than N.
    double normal(int) const;

    /// Get component.
    //  Return the value of the n'th skew multipole component in T/m**(n-1).
    //  Fast version: [b]n[/b] is not checked.
    //  the result is undefined if it is larger than N.
    double skew(int) const;

    /// Get component.
    //  Return a reference to the n'th normal multipole component in T/m**(n-1).
    //  Fast version: [b]n[/b] is not checked.
    //  the result is undefined if it is larger than N.
    double &normal(int);

    /// Get component.
    //  Return a reference to the n'th skew multipole component in T/m**(n-1).
    //  Fast version: [b]n[/b] is not checked.
    //  the result is undefined if it is larger than N.
    double &skew(int);

    /// Add to field.
    //  Add [b]field[/b] to the old field; return new field.
    BMultipoleField &addField(const BMultipoleField &field);

    /// Subtract from field.
    //  Subtract [b]field[/b] from the old field; return new field.
    BMultipoleField &subtractField(const BMultipoleField &field);

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    void scale(double scalar);

    /// Return order.
    int order() const;

private:

    struct Pair {
        // constructors and destructor
        Pair();
        Pair(double normal, double skewed = 0.0);
        Pair(const Pair &);
        ~Pair();

        // operators
        void operator=(const Pair &);
        Pair operator+(const Pair &) const;
        Pair operator-(const Pair &) const;
        Pair operator*(double scale) const;
        void operator+=(const Pair &);
        void operator-=(const Pair &);
        void operator*=(double scale);
        Pair operator-() const;

        double B;   // normal multipole coefficient
        double A;   // skewed multipole coefficient
    };

    // Reserve space for additional coefficients.
    void reserve(int n);

    // The array of multipole coefficients.
    // The dimension is N, where N is the largest order (2N-pole),
    // and pairs[n-1] is the coefficient for the (2n-pole).
    Pair *pairs;

    // The highest order N.
    int itsOrder;
};


// Inline functions
// ------------------------------------------------------------------------

inline int BMultipoleField::order() const {
    return itsOrder;
}


inline double BMultipoleField::getNormalComponent(int n) const {
    if(n > 0 && n <= itsOrder) {
        return pairs[n-1].B;
    } else {
        return 0.0;
    }
}


inline double BMultipoleField::getSkewComponent(int n) const {
    if(n > 0 && n <= itsOrder) {
        return pairs[n-1].A;
    } else {
        return 0.0;
    }
}


inline double &BMultipoleField::normal(int n)
{ return pairs[n-1].B; }


inline double &BMultipoleField::skew(int n)
{ return pairs[n-1].A; }


inline double BMultipoleField::normal(int n) const
{ return pairs[n-1].B; }


inline double BMultipoleField::skew(int n) const
{ return pairs[n-1].A; }

#endif // CLASSIC_BMultipoleField_HH
