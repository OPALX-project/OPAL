// ------------------------------------------------------------------------
// $RCSfile: Octupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: Octupole
//   Defines a concrete representation for a normal (straight) Octupole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/Octupole.h"
#include "Channels/IndirectChannel.h"


// Implementation of typedef Octupole.
// ------------------------------------------------------------------------

// The Octupole type string.
template <>
const std::string SingleMultipole<4>::type("Octupole");

template <>
const SingleMultipole<4>::Entry SingleMultipole<4>::entries[] = {
    {
        "L",
        &SingleMultipole<4>::getElementLength,
        &SingleMultipole<4>::setElementLength
    },
    {
        "B4",
        &SingleMultipole<4>::getComponent,
        &SingleMultipole<4>::setComponent
    },
    { 0, 0, 0 }
};
