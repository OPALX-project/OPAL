//
// Class OpalPepperPot
//   The PEPPERPOT element.
//   The class of OPAL elliptic collimators.
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
#ifndef OPAL_OpalPepperPot_HH
#define OPAL_OpalPepperPot_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalPepperPot: public OpalElement {

public:

    /// The attributes of class OpalPepperPot.
    enum {
        R = COMMON,  // The horizontal half-size of a hole
        NHOLX,
        NHOLY,
        XSIZE,
        YSIZE,
        SIZE
    };

    /// Exemplar constructor.
    OpalPepperPot();

    virtual ~OpalPepperPot();

    /// Make clone.
    virtual OpalPepperPot* clone(const std::string& name);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalPepperPot(const OpalPepperPot&);
    void operator=(const OpalPepperPot&);

    // Clone constructor.
    OpalPepperPot(const std::string& name, OpalPepperPot* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalPepperPot_HH
