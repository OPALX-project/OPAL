#ifndef CLASSIC_OscillatingField_HH
#define CLASSIC_OscillatingField_HH

// ------------------------------------------------------------------------
// $RCSfile: OscillatingField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: OscillatingField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/12/16 19:47:06 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "Fields/EMField.h"
#include <cmath>

// Template class OscillatingField<Field>
// ------------------------------------------------------------------------
/// An oscillating electromagnetic field.
//  The field is derived from the corresponding static field.
//  Care must be taken when the static field is not homogeneous,
//  the oscillating values may be incorrect, as they may not fulfill
//  Maxwell's equations.
//  {p}
//  Note that this class is derived from its template parameter,
//  it inherits the methods and data from that class as well.

template < class Field > class OscillatingField: public Field {

public:

    /// Default constructor.
    //  Constructs a null field.
    OscillatingField();

    virtual ~OscillatingField();

    /// Return the RF frequency in Hz.
    virtual double getFrequency() const;

    /// Return the RF phase in rad.
    virtual double getPhase() const;

    /// Assign the RF frequency in Hz.
    virtual void setFrequency(double f);

    /// Assign the RF phase in rad.
    virtual void setPhase(double phi);

    /// Get field.
    //  Return the constant part of the electric field in point [b]P[/b].
    virtual EVector Efield(const Point3D &point) const;

    /// Get field.
    //  Return the electric field at time [b]t[/b] in point [b]P[/b].
    virtual EVector Efield(const Point3D &point, double time) const;

    /// Get field.
    //  Return the constant part of the magnetic field in point [b]P[/b].
    virtual BVector Bfield(const Point3D &point) const;

    /// Get field.
    //  Return the magnetic field at time [b]t[/b] in point [b]P[/b].
    virtual BVector Bfield(const Point3D &point, double time) const;

private:

    // The frequency, and phase.
    double frequency;
    double phase;
};


// Implementation of template class OscillatingField
// ------------------------------------------------------------------------

template < class Field >
OscillatingField<Field>::OscillatingField(): Field()
{}


template < class Field >
OscillatingField<Field>::~OscillatingField()
{}


template < class Field > inline
double OscillatingField<Field>::getFrequency() const
{ return frequency; }


template < class Field > inline
double OscillatingField<Field>::getPhase() const
{ return phase; }


template < class Field > inline
void OscillatingField<Field>::setFrequency(double freq)
{ frequency = freq; }


template < class Field > inline
void OscillatingField<Field>::setPhase(double phi)
{ phase = phi; }


template < class Field >
EVector OscillatingField<Field>::Efield(const Point3D &p) const
{ return Field::Efield(p); }


template < class Field >
EVector OscillatingField<Field>::Efield(const Point3D &p, double t) const
{ return Field::Efield(p) * cos(t * frequency - phase); }


template < class Field >
BVector OscillatingField<Field>::Bfield(const Point3D &p) const
{ return Field::Bfield(p); }


template < class Field >
BVector OscillatingField<Field>::Bfield(const Point3D &p, double t) const
{ return Field::Bfield(p) * cos(t * frequency - phase); }

#endif // CLASSIC_OscillatingField_HH
