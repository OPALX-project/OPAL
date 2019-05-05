// ------------------------------------------------------------------------
// $RCSfile: CyclotronValleyRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronValleyRep
//   Defines a representation for a Cyclotron Valley.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2010/12/8 14:57:34 $
// $Author: Chuan Wang $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/CyclotronValleyRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(CyclotronValleyRep::*get)() const;
        void (CyclotronValleyRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &CyclotronValleyRep::getElementLength,
            &CyclotronValleyRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class CyclotronValleyRep
// ------------------------------------------------------------------------

CyclotronValleyRep::CyclotronValleyRep():
    CyclotronValley()
{}


CyclotronValleyRep::CyclotronValleyRep(const CyclotronValleyRep &right):
    CyclotronValley(right)
{}


CyclotronValleyRep::CyclotronValleyRep(const std::string &name):
    CyclotronValley(name)
{}


CyclotronValleyRep::~CyclotronValleyRep()
{}


ElementBase *CyclotronValleyRep::clone() const {
    return new CyclotronValleyRep(*this);
}


Channel *CyclotronValleyRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<CyclotronValleyRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}



AcceleratingField &CyclotronValleyRep::getField() {
    return field;
}

const AcceleratingField &CyclotronValleyRep::getField() const {
    return field;
}


StraightGeometry &CyclotronValleyRep::getGeometry() {
    return geometry;
}

const StraightGeometry &CyclotronValleyRep::getGeometry() const {
    return geometry;
}


ElementImage *CyclotronValleyRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}