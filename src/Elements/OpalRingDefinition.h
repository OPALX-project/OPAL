//
// Class OpalRingDefinition
//   The Opal Ring element.
//
// Copyright (c) 2012 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#ifndef OPAL_OpalRingDefinition_HH
#define OPAL_OpalRingDefinition_HH

#include "Elements/OpalElement.h"

class Ring;

/** OpalRingDefinition provides UI wrapper for the Ring
 *
 *  OpalRingDefinition provides User Interface wrapper information for the
 *  Ring. Enables definition of lattice and beam centroid parameters.
 */

class OpalRingDefinition: public OpalElement {
public:
    /** Enumeration maps to UI parameters */
    enum {
        LAT_RINIT = COMMON,
        LAT_PHIINIT,
        LAT_THETAINIT,
        BEAM_RINIT,
        BEAM_PHIINIT,
        BEAM_THETAINIT,
        BEAM_PRINIT,
        HARMONIC_NUMBER,
        SYMMETRY,
        SCALE,
        RFFREQ,
        IS_CLOSED,
        MIN_R,
        MAX_R,
        SIZE // size of the enum
    };

    /** Define mapping from enum variables to string UI parameter names */
    OpalRingDefinition();

    /** No memory allocated so does nothing */
    virtual ~OpalRingDefinition();

    /** Inherited copy constructor */
    virtual OpalRingDefinition* clone(const std::string& name);

    /** Receive parameters from the parser and hand them off to the Ring */
    void update();

    /** Calls print on the OpalElement */
    virtual void print(std::ostream&) const;

private:
    // Not implemented.
    OpalRingDefinition(const OpalRingDefinition&);
    void operator=(const OpalRingDefinition&);

    // Clone constructor.
    OpalRingDefinition(const std::string& name, OpalRingDefinition* parent);
};

#endif // OPAL_OpalRingDefinition_HH
