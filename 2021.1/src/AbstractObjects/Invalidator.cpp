// ------------------------------------------------------------------------
// $RCSfile: Invalidator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Invalidator
//   Abstract base class for references which must be invalidated when an
//   object goes out of scope.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:40:34 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Invalidator.h"


// Class Invalidator
// ------------------------------------------------------------------------


void Invalidator::invalidate() {
    // Default action: do nothing
}
