//
// Class OpalStripper
//   The Stripper element
//
// Copyright (c) 2011, Jianjun Yang, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2014, 2017-2018, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#ifndef OPAL_OpalStripper_HH
#define OPAL_OpalStripper_HH

#include "Elements/OpalElement.h"

class OpalStripper: public OpalElement {

public:

    /// The attributes of class OpalStripper.
    enum {
        XSTART = COMMON,  // Start of x coordinate
        XEND,             // End of x coordinate
        YSTART,           // Start of y coordinate
        YEND,             // End of y coordinate
        WIDTH,            // Width of the probe
        OPCHARGE,         // Charge number of the outcome particle
        OPMASS,           // Mass of the outcome particle
        OPYIELD,
        STOP,
        OUTFN,            // Output filename
        SIZE
    };

    OpalStripper();

    virtual ~OpalStripper();

    /// Make clone.
    virtual OpalStripper *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalStripper(const OpalStripper &);
    void operator=(const OpalStripper &);

    // Clone constructor.
    OpalStripper(const std::string &name, OpalStripper *parent);

};

#endif // OPAL_OpalStripper_HH
