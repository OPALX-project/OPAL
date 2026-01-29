//
// Implements the member function of the EarlyLeaveException class.
//
// Copyright (c) 2008-2019
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#include "Utilities/EarlyLeaveException.h"


EarlyLeaveException::EarlyLeaveException(const std::string &meth, const std::string &msg):
    ClassicException(meth, msg)
{}


EarlyLeaveException::EarlyLeaveException(const EarlyLeaveException &rhs):
    ClassicException(rhs)
{}


EarlyLeaveException::~EarlyLeaveException()
{}
