// ------------------------------------------------------------------------
// $RCSfile: Beamline.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Beamline
//   Represents an abstract sequence of elements.
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/Beamline.h"


// Class Beamline
// ------------------------------------------------------------------------


Beamline::Beamline():
    ElementBase() {
    shareFlag = false;
}


Beamline::Beamline(const Beamline &):
    ElementBase() {
    shareFlag = false;
}


Beamline::Beamline(const std::string &name):
    ElementBase(name) {
    shareFlag = false;
}


Beamline::~Beamline()
{}

Vector_t Beamline::getOrigin3D() const {
    return Vector_t(0);
}

Quaternion Beamline::getInitialDirection() const {
    return Quaternion(1, 0, 0, 0);
}

bool Beamline::getRelativeFlag() const {
    return false;
}