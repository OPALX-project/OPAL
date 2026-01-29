// ------------------------------------------------------------------------
// $RCSfile: Sextupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: Sextupole
//   Defines a concrete representation for a normal (straight) Sextupole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/Sextupole.h"
#include "Channels/IndirectChannel.h"


// Implementation of typedef Sextupole.
// ------------------------------------------------------------------------

// The Sextupole type string.
template <>
const std::string SingleMultipole<3>::type("Sextupole");

template <>
const SingleMultipole<3>::Entry SingleMultipole<3>::entries[] = {
    {
        "L",
        &SingleMultipole<3>::getElementLength,
        &SingleMultipole<3>::setElementLength
    },
    {
        "B3",
        &SingleMultipole<3>::getComponent,
        &SingleMultipole<3>::setComponent
    },
    { 0, 0, 0 }
};
