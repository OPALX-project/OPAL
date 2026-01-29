//
// Class MultipoleRep
//   Representation for a general multipole.
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
#include "BeamlineCore/MultipoleRep.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"
#include <cctype>

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(MultipoleRep::*get)() const;
        void (MultipoleRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &MultipoleRep::getElementLength,
            &MultipoleRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


MultipoleRep::MultipoleRep():
    Multipole(),
    geometry(),
    field()
{}


MultipoleRep::MultipoleRep(const MultipoleRep &multipole):
    Multipole(multipole),
    geometry(multipole.geometry),
    field(multipole.field)
{}


MultipoleRep::MultipoleRep(const std::string &name):
    Multipole(name),
    geometry(),
    field()
{}


MultipoleRep::~MultipoleRep()
{}


ElementBase *MultipoleRep::clone() const {
    return new MultipoleRep(*this);
}


Channel *MultipoleRep::getChannel(const std::string &aKey, bool create) {
    if(aKey[0] == 'A'  ||  aKey[0] == 'B') {
        int n = 0;

        for(std::string::size_type k = 1; k < aKey.length(); k++) {
            if(isdigit(aKey[k])) {
                n = 10 * n + aKey[k] - '0';
            } else {
                return 0;
            }
        }

        if(aKey[0] == 'B') {
            return new IndexedChannel<MultipoleRep>
                   (*this, &MultipoleRep::getNormalComponent,
                    &MultipoleRep::setNormalComponent, n);
        } else {
            return new IndexedChannel<MultipoleRep>
                   (*this, &MultipoleRep::getSkewComponent,
                    &MultipoleRep::setSkewComponent, n);
        }
    } else {
        for(const Entry *entry = entries; entry->name != 0; ++entry) {
            if(aKey == entry->name) {
                return new IndirectChannel<MultipoleRep>
                       (*this, entry->get, entry->set);
            }
        }

        return ElementBase::getChannel(aKey, create);
    }
}


StraightGeometry &MultipoleRep::getGeometry() {
    return geometry;
}

const StraightGeometry &MultipoleRep::getGeometry() const {
    return geometry;
}


BMultipoleField &MultipoleRep::getField() {
    return field;
}


const BMultipoleField &MultipoleRep::getField() const {
    return field;
}


void MultipoleRep::setField(const BMultipoleField &aField) {
    field = aField;
}