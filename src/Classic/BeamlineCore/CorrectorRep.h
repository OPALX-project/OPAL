//
// Class CorrectorRep
//   Representation of a closed orbit corrector.
//   The base class acts on both planes.
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
#ifndef CLASSIC_CorrectorRep_HH
#define CLASSIC_CorrectorRep_HH

#include "AbsBeamline/Corrector.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/BDipoleField.h"

class CorrectorRep: public Corrector {

public:

    /// Constructor with given name.
    explicit CorrectorRep(const std::string &name);

    CorrectorRep();
    CorrectorRep(const CorrectorRep &right);
    virtual ~CorrectorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Get plane(s) of action.
    virtual Plane getPlane() const;

    /// Get horizontal field component in Teslas.
    virtual double getBx() const;

    /// Get vertical field component in Teslas.
    virtual double getBy() const;

    /// Get corrector field.
    virtual BDipoleField &getField();

    /// Get corrector field. Version for const corrector.
    virtual const BDipoleField &getField() const;

    /// Set horizontal field component in Teslas.
    virtual void setBx(double);

    /// Set vertical field component in Teslas.
    virtual void setBy(double);

    /// Set active flag.
    //  If [b]flag[/b] is true, the corrector is activated,
    //  otherwise deactivated.
    virtual void setActive(bool flag = true);

protected:

    /// The corrector geometry.
    StraightGeometry geometry;

    /// The corrector strengths.
    BDipoleField field;

    /// The active/inactive flag.
    bool active;

private:

    // Not implemented.
    void operator=(const CorrectorRep &);
};

#endif // CLASSIC_CorrectorRep_HH
