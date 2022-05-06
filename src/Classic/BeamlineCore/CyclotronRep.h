//
// Class CyclotronRep
//   Representation for a cyclotron magnet system
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_CyclotronRep_HH
#define CLASSIC_CyclotronRep_HH

#include "AbsBeamline/Cyclotron.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"


class CyclotronRep: public Cyclotron {

public:

    /// Constructor with given name.
    explicit CyclotronRep(const std::string &name);

    CyclotronRep();
    CyclotronRep(const CyclotronRep &);
    virtual ~CyclotronRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get the initial radius.
    //  Return the initial radius
    //virtual double getRadius() const;

    /// Set the initial radius
    //  Assign the vertical component of the field in Teslas.
    //  virtual void setRadius(double r);

    /// Get number of slices.
    virtual double getSlices() const;

    /// Get stepsize.
    virtual double getStepsize() const;

    /// Set number of slices.
    virtual void setSlices(double sl);

    /// Set stepsize.
    virtual void setStepsize(double ds);


    /// Set field.
    //  Assign the multipole expansion.
    virtual void setField(const BMultipoleField &field);

    /// Get field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;


    /// Get geometry.
    //  Version for non-constant object.
    virtual PlanarArcGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const PlanarArcGeometry &getGeometry() const;


private:

    // Not implemented.
    void operator=(const CyclotronRep &);

    /// The initial radius of the cyclotron
    double rInit;

    /// The initial momenta of the cyclotron
    double pInit;

    /// The initial phase w.r.t. the rf of the cyclotron
    double phiInit;

    /// The rf of the cyclotron
    double rfFrequ;

    /// The cyclotron geometry.
    PlanarArcGeometry geometry;

    /// The field  expansion.
    BMultipoleField field;

    // Parameters that determine integration step-size.
    double slices;
    double stepsize;
};

#endif // CLASSIC_CyclotronRep_HH
