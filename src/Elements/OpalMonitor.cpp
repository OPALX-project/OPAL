// ------------------------------------------------------------------------
// $RCSfile: OpalMonitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMonitor
//   The class of OPAL monitors for both planes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalMonitor.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MonitorRep.h"
#include "Utilities/Util.h"

// Class OpalMonitor
// ------------------------------------------------------------------------

extern Inform *gmsg;

OpalMonitor::OpalMonitor():
    OpalElement(SIZE, "MONITOR",
                "The \"MONITOR\" element defines a monitor for both planes.") {
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Monitor output filename");

    registerStringAttribute("OUTFN");

    registerOwnership();

    setElement((new MonitorRep("MONITOR"))->makeAlignWrapper());
}


OpalMonitor::OpalMonitor(const std::string &name, OpalMonitor *parent):
    OpalElement(name, parent) {
    setElement((new MonitorRep(name))->makeAlignWrapper());
}


OpalMonitor::~OpalMonitor()
{}


OpalMonitor *OpalMonitor::clone(const std::string &name) {
    return new OpalMonitor(name, this);
}


void OpalMonitor::update() {
    OpalElement::update();

    MonitorRep *mon =
        dynamic_cast<MonitorRep *>(getElement()->removeWrappers());
    double length = std::max(0.01, Attributes::getReal(itsAttr[LENGTH]));
    mon->setElementLength(length);
    mon->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if (Attributes::getString(itsAttr[TYPE]) == "TEMPORAL") {
        mon->setType(Monitor::TEMPORAL);
    } else {
        mon->setType(Monitor::SPATIAL);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mon);
}