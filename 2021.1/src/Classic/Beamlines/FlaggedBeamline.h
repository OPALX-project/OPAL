#ifndef CLASSIC_FlaggedBeamline_HH
#define CLASSIC_FlaggedBeamline_HH

// ------------------------------------------------------------------------
// $RCSfile: FlaggedBeamline.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: FlaggedBeamline
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/FlaggedElmPtr.h"
#include "Beamlines/TBeamline.h"


// Typedef FlaggedBeamline
// ------------------------------------------------------------------------
/// A beam line with flagged elements.
//  An example of a beamline sequence constructed from FlaggedElmPtr
//  objects.  Such a sequence could be used to build a design model with
//  possible flagging of the contained objects.

typedef TBeamline<FlaggedElmPtr> FlaggedBeamline;

#endif // CLASSIC_FlaggedBeamline_HH
