// ------------------------------------------------------------------------
// $RCSfile: SkewQuadrupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SkewQuadrupole
//   Defines a concrete representation for a skewed quadrupole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SkewQuadrupole.h"
#include "Channels/IndirectChannel.h"


// Implementation of typedef SkewQuadrupole.
// ------------------------------------------------------------------------

// The SkewQuadrupole type string.
template <>
const std::string SingleMultipole < -2 >::type("SkewQuadrupole");

template <>
const SingleMultipole < -2 >::Entry SingleMultipole < -2 >::entries[] = {
    {
        "L",
        &SingleMultipole < -2 >::getElementLength,
        &SingleMultipole < -2 >::setElementLength
    },
    {
        "A2",
        &SingleMultipole < -2 >::getComponent,
        &SingleMultipole < -2 >::setComponent
    },
    { 0, 0, 0 }
};
