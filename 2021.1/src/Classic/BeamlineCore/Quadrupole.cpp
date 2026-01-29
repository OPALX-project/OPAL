// ------------------------------------------------------------------------
// $RCSfile: Quadrupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: Quadrupole
//   Defines a concrete representation for a normal (straight) quadrupole.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/Quadrupole.h"
#include "Channels/IndirectChannel.h"


// Implementation of typedef Quadrupole.
// ------------------------------------------------------------------------

// The quadrupole type string.
template <>
const std::string SingleMultipole<2>::type("Quadrupole");

template <>
const SingleMultipole<2>::Entry SingleMultipole<2>::entries[] = {
    {
        "L",
        &SingleMultipole<2>::getElementLength,
        &SingleMultipole<2>::setElementLength
    },
    {
        "B2",
        &SingleMultipole<2>::getComponent,
        &SingleMultipole<2>::setComponent
    },
    { 0, 0, 0 }
};
