//
// This exception class is used to intentionally exit OPAL even if no
// error occurred. It is used in case of the tune calculation with
// the matched distribution option.
//
// Copyright (c) 2008-2019
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#ifndef OPAL_EARLY_LEAVE_EXCEPTION_H
#define OPAL_EARLY_LEAVE_EXCEPTION_H

#include "Utilities/ClassicException.h"

/*!
  This exception class is used to inentionally exit OPAL when no error occurred.
 */
class EarlyLeaveException: public ClassicException {

public:

    /** The usual constructor.
     * @param[in] meth the name of the method or function detecting the exception
     * @param[in] msg the message string identifying the exception
    */
    explicit EarlyLeaveException(const std::string &meth, const std::string &msg);

    EarlyLeaveException(const EarlyLeaveException &);
    virtual ~EarlyLeaveException();

    /// Return the message string for the exception.
    using ClassicException::what;

    /// Return the name of the method or function which detected the exception.
    using ClassicException::where;

private:

    // Not implemented.
    EarlyLeaveException() = delete;
};

#endif
