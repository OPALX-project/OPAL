//
// Class OpalVMonitor
//   The VMONITOR element.
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
#ifndef OPAL_OpalVMonitor_HH
#define OPAL_OpalVMonitor_HH

#include "Elements/OpalElement.h"


class OpalVMonitor: public OpalElement {

public:

    /// Exemplar constructor.
    OpalVMonitor();

    virtual ~OpalVMonitor();

    /// Make clone.
    virtual OpalVMonitor *clone(const std::string &name);

    /// Update the embedded CLASSIC monitor.
    virtual void update();

private:

    // Not implemented.
    OpalVMonitor(const OpalVMonitor &);
    void operator=(const OpalVMonitor &);

    // Clone constructor.
    OpalVMonitor(const std::string &name, OpalVMonitor *parent);
};

#endif // OPAL_OpalVMonitor_HH
