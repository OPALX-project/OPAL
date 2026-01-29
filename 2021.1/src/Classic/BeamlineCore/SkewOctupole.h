#ifndef CLASSIC_SkewOctupole_HH
#define CLASSIC_SkewOctupole_HH

// ------------------------------------------------------------------------
// $RCSfile: SkewOctupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SkewOctupole
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


// Typedef SkewOctupole
// ------------------------------------------------------------------------
/// Representation for a skewed octupole.
//  An instantiation of template SingleMultipole<int>.

typedef SingleMultipole < -4 > SkewOctupole;

#endif // __SkewOctupole_HH
