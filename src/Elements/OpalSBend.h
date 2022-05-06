//
// Class OpalSBend
//   The SBEND element.
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
#ifndef OPAL_OpalSBend_HH
#define OPAL_OpalSBend_HH

#include "Elements/OpalBend.h"

class OpalWake;
class ParticleMatterInteraction;

class OpalSBend: public OpalBend {

public:

    /// Exemplar constructor.
    OpalSBend();

    virtual ~OpalSBend();

    /// Make clone.
    virtual OpalSBend* clone(const std::string& name);

    /// Update the embedded CLASSIC bend.
    virtual void update();

private:

    // Not implemented.
    OpalSBend(const OpalSBend&);
    void operator=(const OpalSBend&);

    // Clone constructor.
    OpalSBend(const std::string& name, OpalSBend* parent);

    OpalWake* owk_m;
    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalSBend_HH