//
// Class StripperRep
//   Representation for Stripper
//
// Copyright (c) 2011, Jianjun Yang,
//                     Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#include "BeamlineCore/StripperRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(StripperRep::*get)() const;
        void (StripperRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &StripperRep::getElementLength,
            &StripperRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class StripperRep
// ------------------------------------------------------------------------

StripperRep::StripperRep():
    Stripper(), field(), geometry(), active(true)
{}


StripperRep::StripperRep(const StripperRep &right):
    Stripper(right), field(), geometry(right.geometry), active(true)
{}


StripperRep::StripperRep(const std::string &name):
    Stripper(name), field(), geometry(), active(true)
{}


StripperRep::~StripperRep()
{}


ElementBase *StripperRep::clone() const {
    return new StripperRep(*this);
}


Channel *StripperRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<StripperRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &StripperRep::getField() {
    return field;
}

const NullField &StripperRep::getField() const {
    return field;
}

StraightGeometry &StripperRep::getGeometry() {
    return geometry;
}

const StraightGeometry &StripperRep::getGeometry() const {
    return geometry;
}


ElementImage *StripperRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}

void StripperRep::setActive(bool flag) {
    active = flag;
}
