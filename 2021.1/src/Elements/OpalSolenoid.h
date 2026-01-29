//
// Class OpalSolenoid
//   The SOLENOID element.
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
#ifndef OPAL_OpalSolenoid_HH
#define OPAL_OpalSolenoid_HH

#include "Elements/OpalElement.h"


class OpalSolenoid: public OpalElement {

public:

    /// The attributes of class OpalSolenoid.
    enum {
        KS = COMMON,  // The longitudinal magnetic field.
        DKS,          // The longitudinal magnetic field error.
        FMAPFN,       // The Field filename
        FAST,         // Faster but less accurate
        SIZE
    };

    /// Exemplar constructor.
    OpalSolenoid();

    virtual ~OpalSolenoid();

    /// Make clone.
    virtual OpalSolenoid *clone(const std::string &name);

    /// Update the embedded CLASSIC solenoid.
    virtual void update();

private:

    // Not implemented.
    OpalSolenoid(const OpalSolenoid &);
    void operator=(const OpalSolenoid &);

    // Clone constructor.
    OpalSolenoid(const std::string &name, OpalSolenoid *parent);
};

#endif // OPAL_OpalSolenoid_HH
