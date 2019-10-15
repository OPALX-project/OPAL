// ------------------------------------------------------------------------
// $RCSfile: BeamStrippingRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamStrippingRep
//   Defines a concrete collimator representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
// ------------------------------------------------------------------------

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