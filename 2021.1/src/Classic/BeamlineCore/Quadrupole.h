#ifndef CLASSIC_Quadrupole_HH
#define CLASSIC_Quadrupole_HH

// ------------------------------------------------------------------------
// $RCSfile: Quadrupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: Quadrupole
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SingleMultipole.h"


// Typedef Quadrupole
// ------------------------------------------------------------------------
/// Representation for a straight (normal) quadrupole.
//  An instantiation of template SingleMultipole<int>.

typedef SingleMultipole<2> Quadrupole;

#endif // __Quadrupole_HH
