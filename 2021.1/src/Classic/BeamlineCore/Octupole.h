#ifndef CLASSIC_Octupole_HH
#define CLASSIC_Octupole_HH

// ------------------------------------------------------------------------
// $RCSfile: Octupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: Octupole
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SingleMultipole.h"


// Typedef Octupole
// ------------------------------------------------------------------------
/// Representation for a straight (normal) Octupole.
//  An instantiation of template SingleMultipole<int>.

typedef SingleMultipole<4> Octupole;

#endif // __Octupole_HH
