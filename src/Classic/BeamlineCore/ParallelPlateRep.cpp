// ------------------------------------------------------------------------
// $RCSfile: ParallelPlateRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelPlateRep
//   Defines a representation for a RF cavity.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/ParallelPlateRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------


namespace {
    struct Entry {
        const char *name;
        double(ParallelPlateRep::*get)() const;
        void (ParallelPlateRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &ParallelPlateRep::getElementLength,
            &ParallelPlateRep::setElementLength
        },
        {
            "AMPLITUDE",
            &ParallelPlateRep::getAmplitude,
            &ParallelPlateRep::setAmplitude
        },
        {
            "FREQUENCY",
            &ParallelPlateRep::getFrequency,
            &ParallelPlateRep::setFrequency
        },
        {
            "PHASE",
            &ParallelPlateRep::getPhase,
            &ParallelPlateRep::setPhase
        },
        { 0, 0, 0 }
    };
}


// Class ParallelPlateRep
// ------------------------------------------------------------------------

ParallelPlateRep::ParallelPlateRep():
    ParallelPlate()
{}


ParallelPlateRep::ParallelPlateRep(const ParallelPlateRep &right):
    ParallelPlate(right)
{}


ParallelPlateRep::ParallelPlateRep(const std::string &name):
    ParallelPlate(name)
{}


ParallelPlateRep::~ParallelPlateRep()
{}


ElementBase *ParallelPlateRep::clone() const {
    return new ParallelPlateRep(*this);
}


Channel *ParallelPlateRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<ParallelPlateRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


AcceleratingField &ParallelPlateRep::getField() {
    return field;
}

const AcceleratingField &ParallelPlateRep::getField() const {
    return field;
}


StraightGeometry &ParallelPlateRep::getGeometry() {
    return geometry;
}

const StraightGeometry &ParallelPlateRep::getGeometry() const {
    return geometry;
}


ElementImage *ParallelPlateRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}