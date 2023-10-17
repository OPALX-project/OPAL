//
// Class OpalMultipoleT
//   The Opal MultipoleT element.
//
// Copyright (c) 2017 - 2023, Titus Dascalu
//                            Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#ifndef OPAL_OPALMULTIPOLET_HH
#define OPAL_OPALMULTIPOLET_HH

#include "Elements/OpalElement.h"

class OpalMultipoleT: public OpalElement {

public:
    // The attributes of class OpalMultipoleT
    enum {
        TP = COMMON,     // Transverse field components
        RFRINGE,         // Length of right fringe field
        LFRINGE,         // Length of left fringe field
        HAPERT,          // Aperture horizontal dimension
        VAPERT,          // Aperture vertical dimension
        MAXFORDER,       // Maximum order in the field expansion
        MAXXORDER,       // Maximum order in x in polynomial expansions
        ANGLE,           // Bending angle of a sector magnet
        ROTATION,        // Rotation angle about central axis for skew elements
        EANGLE,          // Entrance angle
        VARRADIUS,       // Variable radius flag
        BBLENGTH,        // Distance between centre of magnet and entrance
        SIZE             // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalMultipoleT();

    /** Destructor does nothing */
    virtual ~OpalMultipoleT();

    /** Inherited copy constructor */
    virtual OpalMultipoleT* clone(const std::string& name);

    /** Update the MultipoleT with new parameters from UI parser */
    virtual void update();

    void print(std::ostream& os) const;

private:
    // Not implemented.
    OpalMultipoleT(const OpalMultipoleT&);
    void operator=(const OpalMultipoleT&);

    // Clone constructor.
    OpalMultipoleT(const std::string& name, OpalMultipoleT* parent);
};

#endif // OPAL_OpalMultipoleT_HH
