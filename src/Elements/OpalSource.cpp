
#include "Elements/OpalSource.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SourceRep.h"
#include "Physics/Physics.h"


// Class OpalSource
// ------------------------------------------------------------------------

OpalSource::OpalSource():
    OpalElement(SIZE, "SOURCE",
                "The \"SOURCE\" element defines a Source.") {
    itsAttr[DISTRIBUTION] = Attributes::makeStringArray
                             ("DISTRIBUTION", "List of particle distributions to be used ");

    registerStringAttribute("DISTRIBUTION");

    registerOwnership();

    setElement((new SourceRep("SOURCE"))->makeAlignWrapper());
}


OpalSource::OpalSource(const std::string &name, OpalSource *parent):
    OpalElement(name, parent) {
    setElement((new SourceRep(name))->makeAlignWrapper());
}


OpalSource::~OpalSource()
{}


OpalSource *OpalSource::clone(const std::string &name) {
    return new OpalSource(name, this);
}


void OpalSource::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalSource::update() {
    OpalElement::update();

    SourceRep *sol =
        dynamic_cast<SourceRep *>(getElement()->removeWrappers());
    double length = 0.05;

    sol->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sol);
}