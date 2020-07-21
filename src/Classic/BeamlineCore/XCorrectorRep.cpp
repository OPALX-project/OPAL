//
// Class XCorrectorRep
//   Representation for an orbit corrector.
//   This derived class acts on the horizontal plane.
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
#include "BeamlineCore/XCorrectorRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char *name;
        double(XCorrectorRep::*get)() const;
        void (XCorrectorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &XCorrectorRep::getElementLength,
            &XCorrectorRep::setElementLength
        },
        {
            "BY",
            &XCorrectorRep::getBy,
            &XCorrectorRep::setBy
        },
        { 0, 0, 0 }
    };
}


XCorrectorRep::XCorrectorRep():
    CorrectorRep()
{}


XCorrectorRep::XCorrectorRep(const XCorrectorRep &rhs):
    CorrectorRep(rhs)
{}


XCorrectorRep::XCorrectorRep(const std::string &name):
    CorrectorRep(name)
{}


XCorrectorRep::~XCorrectorRep()
{}


ElementBase *XCorrectorRep::clone() const {
    return new XCorrectorRep(*this);
}


Channel *XCorrectorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<XCorrectorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Corrector::Plane XCorrectorRep::getPlane() const {
    return active ? X : OFF;
}


double XCorrectorRep::getBx() const {
    return 0.0;
}


void XCorrectorRep::setBx(double) {
    // Do nothing.
}
