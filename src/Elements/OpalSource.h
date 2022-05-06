//
// Class OpalSource
//   The SOURCE element.
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
#ifndef OPAL_OPALSOURCE_HH
#define OPAL_OPALSOURCE_HH

#include "Elements/OpalElement.h"


class OpalSource: public OpalElement {

public:

    /// The attributes of class OpalSource.
    enum {
        DISTRIBUTION = COMMON,  // The longitudinal magnetic field.
        TRANSPARENT,
        SIZE
    };

    /// Exemplar constructor.
    OpalSource();

    virtual ~OpalSource();

    /// Make clone.
    virtual OpalSource* clone(const std::string& name);

    /// Update the embedded CLASSIC solenoid.
    virtual void update();

private:

    // Not implemented.
    OpalSource(const OpalSource&);
    void operator=(const OpalSource&);

    // Clone constructor.
    OpalSource(const std::string& name, OpalSource* parent);
};

#endif // OPAL_OPALSOURCE_HH
