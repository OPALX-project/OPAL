//
// Class XCorrectorRep
//   Representation for an orbit corrector.
//   This derived class acts on the horizontal plane.
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
#ifndef CLASSIC_XCorrectorRep_HH
#define CLASSIC_XCorrectorRep_HH

#include "BeamlineCore/CorrectorRep.h"


class XCorrectorRep: public CorrectorRep {

public:

    /// Constructor with given name.
    explicit XCorrectorRep(const std::string &name);

    XCorrectorRep();
    XCorrectorRep(const XCorrectorRep &);
    virtual ~XCorrectorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get plane of action.
    //  Return the x-plane for this class.
    virtual Plane getPlane() const;

    /// Get field.
    //  Return horizontal component (always zero).
    virtual double getBx() const;

    /// Set field.
    //  Ignore the horizontal field value.
    virtual void setBx(double);

private:

    // Not implemented.
    void operator=(const XCorrectorRep &);
};

#endif // CLASSIC_XCorrectorRep_HH
