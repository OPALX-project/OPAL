#ifndef CLASSIC_SkewSextupole_HH
#define CLASSIC_SkewSextupole_HH

// ------------------------------------------------------------------------
// $RCSfile: SkewSextupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SkewSextupole
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


// Typedef SkewSextupole
// ------------------------------------------------------------------------
/// Representation for a skewed sextupole.
//  An instantiation of template SingleMultipole<int>.

typedef SingleMultipole < -3 > SkewSextupole;

#endif // __SkewSextupole_HH
