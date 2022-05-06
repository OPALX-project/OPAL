//
// Class OpalMonitor
//   The MONITOR element.
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
#ifndef OPAL_OpalMonitor_HH
#define OPAL_OpalMonitor_HH

#include "Elements/OpalElement.h"

class OpalMonitor: public OpalElement {

public:

    /// Exemplar constructor.
    OpalMonitor();

    virtual ~OpalMonitor();

    /// Make clone.
    virtual OpalMonitor* clone(const std::string& name);

    /// Update the embedded CLASSIC monitor.
    virtual void update();

private:

    // Not implemented.
    OpalMonitor(const OpalMonitor&);
    void operator=(const OpalMonitor&);

    // Clone constructor.
    OpalMonitor(const std::string& name, OpalMonitor* parent);
};

#endif // OPAL_OpalMonitor_HH
