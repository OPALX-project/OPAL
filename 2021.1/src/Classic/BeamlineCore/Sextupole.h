#ifndef CLASSIC_Sextupole_HH
#define CLASSIC_Sextupole_HH

// ------------------------------------------------------------------------
// $RCSfile: Sextupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: Sextupole
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


// Typedef Sextupole
// ------------------------------------------------------------------------
/// Representation for a straight (normal) Sextupole.
//  An instantiation of template SingleMultipole<int>.

typedef SingleMultipole<3> Sextupole;

#endif // __Sextupole_HH
