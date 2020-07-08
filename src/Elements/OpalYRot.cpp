// ------------------------------------------------------------------------
// $RCSfile: OpalYRot.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalYRot
//   The class of OPAL YROT coordinate transforms.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalYRot.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/PatchRep.h"


// Class OpalYRot
// ------------------------------------------------------------------------

OpalYRot::OpalYRot():
    OpalElement(SIZE, "YROT",
                "The \"YROT\" element defines a coordinate rotation "
                "around the y-axis.") {
    itsAttr[ANGLE] = Attributes::makeReal
                     ("ANGLE", "Angle for the transformation");

    registerRealAttribute("YANGLE");

    registerOwnership();

    setElement(new PatchRep("YROT"));
}


OpalYRot::OpalYRot(const std::string &name, OpalYRot *parent):
    OpalElement(name, parent) {
    setElement(new PatchRep(name));
}


OpalYRot::~OpalYRot()
{}


OpalYRot *OpalYRot::clone(const std::string &name) {
    return new OpalYRot(name, this);
}


void OpalYRot::
fillRegisteredAttributes(const ElementBase &base) {
    OpalElement::fillRegisteredAttributes(base);
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    attributeRegistry["SANGLE"]->setReal(angle);
}


void OpalYRot::update() {
    // Define geometry.
    PatchRep *yrot = dynamic_cast<PatchRep *>(getElement());
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    Euclid3D rot = Euclid3D::YRotation(- angle);
    yrot->setPatch(rot);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(yrot);
}