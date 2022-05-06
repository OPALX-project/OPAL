//
// Class OpalSlit
//   The SLIT element.
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
#ifndef OPAL_OpalSlit_HH
#define OPAL_OpalSlit_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalSlit: public OpalElement {

public:

    /// The attributes of class OpalSlit.
    enum {
        XSIZE = COMMON,  // The horizontal half-size.
        YSIZE,           // The vertical half-size.
        SIZE
    };

    /// Exemplar constructor.
    OpalSlit();

    virtual ~OpalSlit();

    /// Make clone.
    virtual OpalSlit* clone(const std::string& name);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalSlit(const OpalSlit&);
    void operator=(const OpalSlit&);

    // Clone constructor.
    OpalSlit(const std::string& name, OpalSlit* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalSlit_HH
