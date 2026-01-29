// ------------------------------------------------------------------------
// $RCSfile: SkewOctupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SkewOctupole
//   Defines a concrete representation for a skewed octupole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SkewOctupole.h"
#include "Channels/IndirectChannel.h"


// Implementation of typedef SkewOctupole.
// ------------------------------------------------------------------------

// The SkewOctupole type string.
template <>
const std::string SingleMultipole < -4 >::type("SkewOctupole");

template <>
const SingleMultipole < -4 >::Entry SingleMultipole < -4 >::entries[] = {
    {
        "L",
        &SingleMultipole < -4 >::getElementLength,
        &SingleMultipole < -4 >::setElementLength
    },
    {
        "A4",
        &SingleMultipole < -4 >::getComponent,
        &SingleMultipole < -4 >::setComponent
    },
    { 0, 0, 0 }
};
