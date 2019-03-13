// ------------------------------------------------------------------------
// $RCSfile: UndulatorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: UndulatorRep
//   Defines a concrete drift space representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/UndulatorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(UndulatorRep::*get)() const;
        void (UndulatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &UndulatorRep::getElementLength,
            &UndulatorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class UndulatorRep
// ------------------------------------------------------------------------

UndulatorRep::UndulatorRep():
    Undulator(),
    geometry(0.0),
    k_m(0.0),
    fname_m("")
{}


UndulatorRep::UndulatorRep(const UndulatorRep &right):
    Undulator(right),
    geometry(right.geometry),
    k_m(right.k_m),
    fname_m(right.fname_m)
{}


UndulatorRep::UndulatorRep(const std::string &name):
    Undulator(name),
    geometry(),
    k_m(0.0),
    fname_m("")
{}


UndulatorRep::~UndulatorRep()
{}


ElementBase *UndulatorRep::clone() const {
    return new UndulatorRep(*this);
}


Channel *UndulatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<UndulatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &UndulatorRep::getField() {
    return field;
}

const NullField &UndulatorRep::getField() const {
    return field;
}


StraightGeometry &UndulatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &UndulatorRep::getGeometry() const {
    return geometry;
}


ElementImage *UndulatorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}

/*
  Here we 

 */

void   UndulatorRep::setK(double k) { k_m = k; }
double UndulatorRep::getK() { return k_m; }

void UndulatorRep::setFilename(const std::string&fname) {
  fname_m = fname;
}

const std::string& UndulatorRep::getFilename() const {
  return fname_m;
}
