//
// Class YCorrectorRep
//   Representation for an orbit corrector.
//   Acts on the vertical plane.
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
#include "BeamlineCore/YCorrectorRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char *name;
        double(YCorrectorRep::*get)() const;
        void (YCorrectorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &YCorrectorRep::getElementLength,
            &YCorrectorRep::setElementLength
        },
        {
            "BX",
            &YCorrectorRep::getBx,
            &YCorrectorRep::setBx
        },
        { 0, 0, 0 }
    };
}


YCorrectorRep::YCorrectorRep():
    CorrectorRep()
{}


YCorrectorRep::YCorrectorRep(const YCorrectorRep &right):
    CorrectorRep(right)
{}


YCorrectorRep::YCorrectorRep(const std::string &name):
    CorrectorRep(name)
{}


YCorrectorRep::~YCorrectorRep()
{}


ElementBase *YCorrectorRep::clone() const {
    return new YCorrectorRep(*this);
}


Channel *YCorrectorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<YCorrectorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


Corrector::Plane YCorrectorRep::getPlane() const {
    return Y;
}


double YCorrectorRep::getBy() const {
    return 0.0;
}


void YCorrectorRep::setBy(double) {
    // Do nothing.
}
