//
// Class OpalFlexibleCollimator
//   The Flexible Collimator element.
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
#include "Elements/OpalFlexibleCollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Utilities/OpalException.h"

#include <boost/regex.hpp>

OpalFlexibleCollimator::OpalFlexibleCollimator():
    OpalElement(SIZE, "FLEXIBLECOLLIMATOR",
                "The \"FLEXIBLECOLLIMATOR\" element defines a flexible collimator."),
    parmatint_m(nullptr) {
    itsAttr[FNAME] = Attributes::makeString
                     ("FNAME", "File name containing description of holes");
    itsAttr[DESC]  = Attributes::makeString
                     ("DESCRIPTION", "String describing the distribution of holes");
    itsAttr[DUMP]  = Attributes::makeBool
                     ("DUMP", "Save quadtree and holes of collimator", false);
    registerOwnership();

    setElement(new FlexibleCollimatorRep("FLEXIBLECOLLIMATOR"));
}


OpalFlexibleCollimator::OpalFlexibleCollimator(const std::string& name, OpalFlexibleCollimator* parent):
    OpalElement(name, parent),
    parmatint_m(nullptr) {
    setElement(new FlexibleCollimatorRep(name));
}


OpalFlexibleCollimator::~OpalFlexibleCollimator() {
    if (parmatint_m)
        delete parmatint_m;
}


OpalFlexibleCollimator* OpalFlexibleCollimator::clone(const std::string& name) {
    return new OpalFlexibleCollimator(name, this);
}


void OpalFlexibleCollimator::update() {
    OpalElement::update();

    FlexibleCollimatorRep* coll =
        dynamic_cast<FlexibleCollimatorRep*>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);

    std::string fname = Attributes::getString(itsAttr[FNAME]);
    std::string desc = Attributes::getString(itsAttr[DESC]);
    if (!fname.empty()) {
        std::ifstream it(fname);
        std::string str((std::istreambuf_iterator<char>(it)),
                        std::istreambuf_iterator<char>());

        str = boost::regex_replace(str, boost::regex("//.*?\\n"), std::string(""), boost::match_default | boost::format_all);
        str = boost::regex_replace(str, boost::regex("\\s"), std::string(""), boost::match_default | boost::format_all);

        coll->setDescription(str);
    } else if (!desc.empty()) {
        desc = boost::regex_replace(desc, boost::regex("[\\t ]"), std::string(""), boost::match_default | boost::format_all);
        coll->setDescription(desc);
    } else if (getOpalName() != "FLEXIBLECOLLIMATOR") {
        throw OpalException("OpalFlexibleCollimator::update",
                            "A description for the holes has to be provided, either using DESCRIPTION or FNAME");
    }
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if (itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == nullptr) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    if (Attributes::getBool(itsAttr[DUMP])) {
        coll->writeHolesAndQuadtree(getOpalName());
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}
