// ------------------------------------------------------------------------
// $RCSfile: ConstChannel.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstChannel
//   A channel which gives read only access to a variable.
//
// ------------------------------------------------------------------------
// Class category: Channels
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Channels/ConstChannel.h"


// Class ConstChannel
// ------------------------------------------------------------------------


ConstChannel::ConstChannel()
{}


ConstChannel::~ConstChannel()
{}


bool ConstChannel::isSettable() const {
    // Default: ConstChannel is never settable.
    return false;
}
