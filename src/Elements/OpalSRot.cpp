// ------------------------------------------------------------------------
// $RCSfile: OpalSRot.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSRot
//   The class of OPAL SROT coordinate transforms.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSRot.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/PatchRep.h"


// Class OpalSRot
// ------------------------------------------------------------------------

OpalSRot::OpalSRot():
    OpalElement(SIZE, "SROT",
                "The \"SRot\" element defines a coordinate rotation "
                "around the s-axis.") {
    itsAttr[ANGLE] = Attributes::makeReal
                     ("ANGLE", "Angle for the transformation");

    registerRealAttribute("ANGLE");

    registerOwnership();

    setElement(new PatchRep("SROT"));
}


OpalSRot::OpalSRot(const std::string &name, OpalSRot *parent):
    OpalElement(name, parent) {
    setElement(new PatchRep(name));
}


OpalSRot::~OpalSRot()
{}


OpalSRot *OpalSRot::clone(const std::string &name) {
    return new OpalSRot(name, this);
}


void OpalSRot::
fillRegisteredAttributes(const ElementBase &base) {
    OpalElement::fillRegisteredAttributes(base);

    double angle = Attributes::getReal(itsAttr[ANGLE]);
    attributeRegistry["SANGLE"]->setReal(angle);
}


void OpalSRot::update() {
    // Define geometry.
    PatchRep *srot = dynamic_cast<PatchRep *>(getElement());
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    Euclid3D rot = Euclid3D::ZRotation(angle);
    srot->setPatch(rot);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(srot);
}