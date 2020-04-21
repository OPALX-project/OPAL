//
// Class BeamStrippingRep
//   Defines a concrete representation for beam stripping.
//
// Copyright (c) 2018-2019, Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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

#include "AbsBeamline/ElementImage.h"
#include "BeamlineCore/BeamStrippingRep.h"
#include "Channels/IndirectChannel.h"

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(BeamStrippingRep::*get)() const;
        void (BeamStrippingRep::*set)(double);
    };
    const Entry entries[] = {
        {
            "p",
            &BeamStripping::getPressure,
            &BeamStripping::setPressure
        },
        { 0, 0, 0 }
    };
}


// Class BeamStrippingRep
// ------------------------------------------------------------------------

BeamStrippingRep::BeamStrippingRep():
    BeamStripping(),
    geometry(0.0, 0.0)
{}

BeamStrippingRep::BeamStrippingRep(const BeamStrippingRep &right):
    BeamStripping(right),
    geometry(0.0, 0.0)
{}

BeamStrippingRep::BeamStrippingRep(const std::string &name):
    BeamStripping(name),
    geometry(0.0, 0.0)
{}

BeamStrippingRep::~BeamStrippingRep()
{}


ElementBase *BeamStrippingRep::clone() const {
    return new BeamStrippingRep(*this);
}


Channel *BeamStrippingRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<BeamStrippingRep>(*this, entry->get, entry->set);
        }
    }
    return ElementBase::getChannel(aKey, create);
}


NullField &BeamStrippingRep::getField() {
    return field;
}

const NullField &BeamStrippingRep::getField() const {
    return field;
}


PlanarArcGeometry &BeamStrippingRep::getGeometry() {
    return geometry;
}

const PlanarArcGeometry &BeamStrippingRep::getGeometry() const {
    return geometry;
}


ElementImage *BeamStrippingRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}