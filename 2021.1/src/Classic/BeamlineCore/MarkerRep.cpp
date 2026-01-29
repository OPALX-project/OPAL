// ------------------------------------------------------------------------
// $RCSfile: MarkerRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MarkerRep
//   Defines a concrete marker element.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/MarkerRep.h"


// Class MarkerRep
// ------------------------------------------------------------------------

MarkerRep::MarkerRep():
    Marker()
{}

MarkerRep::MarkerRep(const MarkerRep &right):
    Marker(right)
{}


MarkerRep::MarkerRep(const std::string &name):
    Marker(name)
{}


MarkerRep::~MarkerRep()
{}


ElementBase *MarkerRep::clone() const {
    return new MarkerRep(*this);
}


NullField &MarkerRep::getField() {
    return field;
}

const NullField &MarkerRep::getField() const {
    return field;
}


NullGeometry &MarkerRep::getGeometry() {
    return geometry;
}

const NullGeometry &MarkerRep::getGeometry() const {
    return geometry;
}


double MarkerRep::getArcLength() const {
    return 0.0;
}


double MarkerRep::getElementLength() const {
    return 0.0;
}
