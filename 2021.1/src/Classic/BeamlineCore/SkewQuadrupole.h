#ifndef CLASSIC_SkewQuadrupole_HH
#define CLASSIC_SkewQuadrupole_HH

// ------------------------------------------------------------------------
// $RCSfile: SkewQuadrupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SkewQuadrupole
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


// Typedef SkewQuadrupole
// ------------------------------------------------------------------------
/// Representation for a skewed quadrupole.
//  An instantiation of template SingleMultipole<int>.

typedef SingleMultipole < -2 > SkewQuadrupole;

#endif // __SkewQuadrupole_HH
