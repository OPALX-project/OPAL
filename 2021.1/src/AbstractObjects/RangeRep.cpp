// ------------------------------------------------------------------------
// $RCSfile: RangeRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RangeRep
//   A class used to represent a range specification.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:05:47 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/RangeRep.h"
#include "AbstractObjects/Element.h"
#include <iostream>


// Class RangeRep
// ------------------------------------------------------------------------

RangeRep::RangeRep():
    first(), last(), fullRange(true)
{}


RangeRep::RangeRep(const RangeRep &rep):
    first(rep.first), last(rep.last), fullRange(rep.fullRange)
{}


RangeRep::RangeRep(PlaceRep &fst, PlaceRep &lst):
    first(fst), last(lst), fullRange(false)
{}


RangeRep::~RangeRep()
{}


const RangeRep &RangeRep::operator=(const RangeRep &rep) {
    first     = rep.first;
    last      = rep.last;
    fullRange = rep.fullRange;
    return *this;
}


void RangeRep::initialize() {
    if(! fullRange) {
        first.initialize();
        last.initialize();
        status = false;
    } else {
        status = true;
    }
}


bool RangeRep::isActive() const {
    return status;
}


void RangeRep::enter(const FlaggedElmPtr &fep) const {
    if(! fullRange) {
        // Enter range, if we are in first element.
        first.enter(fep);
        if(first.isActive()) status = true;
        last.enter(fep);
    }
}


void RangeRep::leave(const FlaggedElmPtr &fep) const {
    if(! fullRange) {
        // Leave range, if we are in last element.
        if(last.isActive()) status = false;
        first.leave(fep);
        last.leave(fep);
    }
}


void RangeRep::print(std::ostream &os) const {
    if(fullRange) {
        os << "FULL";
    } else {
        first.print(os);
        os << '/';
        last.print(os);
    }

    return;
}
