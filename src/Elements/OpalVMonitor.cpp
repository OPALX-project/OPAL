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
#include "Elements/OpalVMonitor.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/YMonitorRep.h"


OpalVMonitor::OpalVMonitor():
    OpalElement(COMMON, "VMONITOR",
                "The \"VMONITOR\" element defines a monitor "
                "for the vertical plane.") {
    setElement(new YMonitorRep("VMONITOR"));
}


OpalVMonitor::OpalVMonitor(const std::string &name, OpalVMonitor *parent):
    OpalElement(name, parent) {
    setElement(new YMonitorRep(name));
}


OpalVMonitor::~OpalVMonitor()
{}


OpalVMonitor *OpalVMonitor::clone(const std::string &name) {
    return new OpalVMonitor(name, this);
}


void OpalVMonitor::update() {
    OpalElement::update();

    YMonitorRep *mon =
        dynamic_cast<YMonitorRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    mon->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}
