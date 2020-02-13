// ------------------------------------------------------------------------
// $RCSfile: ConstBField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstBField
//   A static magnetic field independent of (x,y,z).
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstBField.h"


// Class ConstBField
// ------------------------------------------------------------------------

ConstBField::ConstBField()
{}


ConstBField::~ConstBField()
{}


double ConstBField::getBx() const {
    return 0.0;
}

double ConstBField::getBy() const {
    return 0.0;
}

double ConstBField::getBz() const {
    return 0.0;
}


void ConstBField::setBx(double /*B*/)
{}


void ConstBField::setBy(double /*B*/)
{}


void ConstBField::setBz(double /*B*/)
{}
