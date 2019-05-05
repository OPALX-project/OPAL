// ------------------------------------------------------------------------
// $RCSfile: DegraderRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DegraderRep
//   Defines a concrete collimator representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/DegraderRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(DegraderRep::*get)() const;
        void (DegraderRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &DegraderRep::getElementLength,
            &DegraderRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class DegraderRep
// ------------------------------------------------------------------------

DegraderRep::DegraderRep():
    Degrader(),
    geometry(0.0)
{}


DegraderRep::DegraderRep(const DegraderRep &right):
    Degrader(right),
    geometry(right.geometry)
{}


DegraderRep::DegraderRep(const std::string &name):
    Degrader(name),
    geometry()
{}


DegraderRep::~DegraderRep()
{}


ElementBase *DegraderRep::clone() const {
    return new DegraderRep(*this);
}


Channel *DegraderRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<DegraderRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &DegraderRep::getField() {
    return field;
}

const NullField &DegraderRep::getField() const {
    return field;
}


StraightGeometry &DegraderRep::getGeometry() {
    return geometry;
}

const StraightGeometry &DegraderRep::getGeometry() const {
    return geometry;
}


ElementImage *DegraderRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}