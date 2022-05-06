//
// Class YCorrectorRep
//   Representation for an orbit corrector.
//   Acts on the vertical plane.
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
#ifndef CLASSIC_YCorrectorRep_HH
#define CLASSIC_YCorrectorRep_HH

#include "BeamlineCore/CorrectorRep.h"


class YCorrectorRep: public CorrectorRep {

public:

    /// Constructor with given name.
    explicit YCorrectorRep(const std::string &name);

    YCorrectorRep();
    YCorrectorRep(const YCorrectorRep &);
    virtual ~YCorrectorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get plane.
    //  Return the y-plane for this class.
    virtual Plane getPlane() const;

    /// Get field.
    //  Return the vertical field component (always zero).
    virtual double getBy() const;

    /// Set field.
    //  Ignore the vertical field component.
    virtual void setBy(double);

private:

    // Not implemented.
    void operator=(const YCorrectorRep &);
};

#endif // CLASSIC_YCorrectorRep_HH
