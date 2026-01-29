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
#include "Elements/OpalMonitor.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MonitorRep.h"
#include "Utilities/Util.h"


OpalMonitor::OpalMonitor():
    OpalElement(COMMON, "MONITOR",
                "The \"MONITOR\" element defines a monitor for both planes.") {

    registerOwnership();

    setElement(new MonitorRep("MONITOR"));
}


OpalMonitor::OpalMonitor(const std::string& name, OpalMonitor* parent):
    OpalElement(name, parent) {
    setElement(new MonitorRep(name));
}


OpalMonitor::~OpalMonitor()
{}


OpalMonitor* OpalMonitor::clone(const std::string& name) {
    return new OpalMonitor(name, this);
}


void OpalMonitor::update() {
    OpalElement::update();

    MonitorRep* mon =
        dynamic_cast<MonitorRep*>(getElement());

    double length = std::max(0.01, Attributes::getReal(itsAttr[LENGTH]));
    mon->setElementLength(length);
    mon->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if (Attributes::getString(itsAttr[TYPE]) == "TEMPORAL") {
        mon->setCollectionType(CollectionType::TEMPORAL);
    } else {
        mon->setCollectionType(CollectionType::SPATIAL);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}
