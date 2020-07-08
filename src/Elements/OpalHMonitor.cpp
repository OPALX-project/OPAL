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
#include "Elements/OpalHMonitor.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/XMonitorRep.h"


OpalHMonitor::OpalHMonitor():
    OpalElement(COMMON, "HMONITOR",
                "The \"HMONITOR\" element defines a monitor "
                "for the horizontal plane.") {
    setElement(new XMonitorRep("HMONITOR"));
}


OpalHMonitor::OpalHMonitor(const std::string &name, OpalHMonitor *parent):
    OpalElement(name, parent) {
    setElement(new XMonitorRep(name));
}


OpalHMonitor::~OpalHMonitor()
{}


OpalHMonitor *OpalHMonitor::clone(const std::string &name) {
    return new OpalHMonitor(name, this);
}


void OpalHMonitor::update() {
    OpalElement::update();

    XMonitorRep *mon =
        dynamic_cast<XMonitorRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    mon->setElementLength(length);

   // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}
