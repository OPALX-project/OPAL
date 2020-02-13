// ------------------------------------------------------------------------
// $RCSfile: ConstEField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstEField
//   A static electric field independent of (x,y,z).
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstEField.h"


// Class ConstEField
// ------------------------------------------------------------------------

ConstEField::ConstEField()
{}


ConstEField::~ConstEField()
{}


double ConstEField::getEx() const {
    return 0.0;
}


double ConstEField::getEy() const {
    return 0.0;
}


double ConstEField::getEz() const {
    return 0.0;
}


void ConstEField::setEx(double /*E*/)
{}


void ConstEField::setEy(double /*E*/)
{}


void ConstEField::setEz(double /*E*/)
{}
