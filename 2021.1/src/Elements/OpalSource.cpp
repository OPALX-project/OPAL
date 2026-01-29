//
// Class OpalSource
//   The SOURCE element.
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
#include "Elements/OpalSource.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SourceRep.h"
#include "Physics/Physics.h"


OpalSource::OpalSource():
    OpalElement(SIZE, "SOURCE",
                "The \"SOURCE\" element defines a Source.") {
    itsAttr[DISTRIBUTION] = Attributes::makeStringArray
        ("DISTRIBUTION", "List of particle distributions to be used ");

    itsAttr[TRANSPARENT] = Attributes::makeBool
        ("TRANSPARENT", "Make the source element transparent to impacting elements; Default value is FALSE", false);

    registerOwnership();

    setElement(new SourceRep("SOURCE"));
}


OpalSource::OpalSource(const std::string& name, OpalSource* parent):
    OpalElement(name, parent) {
    setElement(new SourceRep(name));
}


OpalSource::~OpalSource()
{}


OpalSource *OpalSource::clone(const std::string& name) {
    return new OpalSource(name, this);
}


void OpalSource::update() {
    OpalElement::update();

    SourceRep* sol =
        dynamic_cast<SourceRep*>(getElement());
    double length = 0.05;

    sol->setElementLength(length);

    if (Attributes::getBool(itsAttr[TRANSPARENT])) {
        sol->setTransparent();
    }

    sol->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sol);
}
