// ------------------------------------------------------------------------
// $RCSfile: NullField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NullField
//   A class for zero electromagnetic fields.
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:36 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/NullField.h"


// Class NullField
// ------------------------------------------------------------------------

NullField::NullField()
{}


NullField::NullField(const NullField & rhs):
  EMField(rhs)
{}


NullField::~NullField()
{}


const NullField &NullField::operator=(const NullField &)
{ return *this; }


void NullField::scale(double)
{}
