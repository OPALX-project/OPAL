// ------------------------------------------------------------------------
// $RCSfile: CyclotronRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronRep
//   Defines a concrete representation for a cyclotron.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/CyclotronRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"
#include "ComponentWrappers/CyclotronWrapper.h"

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(CyclotronRep::*get)() const;
        void (CyclotronRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "RINIT",
            0, // :FIXME: Why commented out? &CyclotronRep::getRadius, 
            0  // :FIXME: &CyclotronRep::setRadius
        },
        { 0, 0, 0 }
    };
}


// Class CyclotronRep
// ------------------------------------------------------------------------

CyclotronRep::CyclotronRep():
    Cyclotron(),
    geometry(0.0, 0.0),
    field() {
    rInit = 0.0;
}


CyclotronRep::CyclotronRep(const CyclotronRep &rhs):
    Cyclotron(rhs),
    geometry(rhs.geometry),
    field(rhs.field) {
    rInit = rhs.rInit;
}


CyclotronRep::CyclotronRep(const std::string &name):
    Cyclotron(name),
    geometry(0.0, 0.0),
    field() {
    rInit = 0.0;
}


CyclotronRep::~CyclotronRep()
{}


ElementBase *CyclotronRep::clone() const {
    return new CyclotronRep(*this);
}


Channel *CyclotronRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
      if(aKey == entry->name) {
        return new IndirectChannel<CyclotronRep>(*this, entry->get, entry->set);
      }
    }

    return ElementBase::getChannel(aKey, create);
}


double CyclotronRep::getSlices() const {
    return slices;
}

double CyclotronRep::getStepsize() const {
    return stepsize;
}

void CyclotronRep::setSlices(double sl) {
    slices = sl;
}

void CyclotronRep::setStepsize(double ds) {
    stepsize = ds;
}

/*
void CyclotronRep::setRadius(double r)
{
  rInit = r;
}

double CyclotronRep::getRadius() const
{
  return rInit ;
}
*/

PlanarArcGeometry &CyclotronRep::getGeometry() {
    return geometry;
}

const PlanarArcGeometry &CyclotronRep::getGeometry() const {
    return geometry;
}

BMultipoleField &CyclotronRep::getField() {
    return field;
}

const BMultipoleField &CyclotronRep::getField() const {
    return field;
}

void CyclotronRep::setField(const BMultipoleField &f) {
    field = f;
}

ElementBase *CyclotronRep::makeFieldWrapper() {
    ElementBase *wrap = new CyclotronWrapper(this);
    wrap->setName(getName());
    return wrap;
}

