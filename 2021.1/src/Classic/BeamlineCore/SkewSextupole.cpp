// ------------------------------------------------------------------------
// $RCSfile: SkewSextupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SkewSextupole
//   Defines a concrete representation for a skewed sextupole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SkewSextupole.h"
#include "Channels/IndirectChannel.h"


// Implementation of typedef SkewSextupole.
// ------------------------------------------------------------------------

// The SkewSextupole type string.
template <>
const std::string SingleMultipole < -3 >::type("SkewSextupole");

template <>
const SingleMultipole < -3 >::Entry SingleMultipole < -3 >::entries[] = {
    {
        "L",
        &SingleMultipole < -3 >::getElementLength,
        &SingleMultipole < -3 >::setElementLength
    },
    {
        "A3",
        &SingleMultipole < -3 >::getComponent,
        &SingleMultipole < -3 >::setComponent
    },
    { 0, 0, 0 }
};
