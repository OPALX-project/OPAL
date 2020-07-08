//
// Class OpalHMonitor
//   The HMONITOR element.
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
#ifndef OPAL_OpalHMonitor_HH
#define OPAL_OpalHMonitor_HH

#include "Elements/OpalElement.h"

class OpalHMonitor: public OpalElement {

public:

    /// Exemplar constructor.
    OpalHMonitor();

    virtual ~OpalHMonitor();

    /// Make clone.
    virtual OpalHMonitor *clone(const std::string &name);

    /// Update the embedded CLASSIC monitor.
    virtual void update();

private:

    // Not implemented.
    OpalHMonitor(const OpalHMonitor &);
    void operator=(const OpalHMonitor &);

    // Clone constructor.
    OpalHMonitor(const std::string &name, OpalHMonitor *parent);
};

#endif // OPAL_OpalHMonitor_HH
