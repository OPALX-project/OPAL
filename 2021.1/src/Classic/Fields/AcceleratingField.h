#ifndef CLASSIC_AcceleratingField_HH
#define CLASSIC_AcceleratingField_HH

// ------------------------------------------------------------------------
// $RCSfile: AcceleratingField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AcceleratingField (typedef)
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/ConstEzField.h"
#include "Fields/OscillatingField.h"


// Typedef AcceleratingField
// ------------------------------------------------------------------------
/// The electromagnetic field of an RF cavity.
//  A homogeneous field in z-direction, oscillating in time.
//  An instance of the template class OscillatingField.

typedef OscillatingField<ConstEzField> AcceleratingField;

#endif // CLASSIC_AcceleratingField_HH
