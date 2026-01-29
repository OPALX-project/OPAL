#ifndef CLASSIC_ConvergenceError_HH
#define CLASSIC_ConvergenceError_HH

// ------------------------------------------------------------------------
// $RCSfile: ConvergenceError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConvergenceError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ArithmeticError.h"

#include <string>

// Class ConvergenceError
// ------------------------------------------------------------------------
/// Convergence error exception.
//  This exception is thrown, when an iterative algorithm fails to converge.

class ConvergenceError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    ConvergenceError(const std::string &method, const std::string &message);

    ConvergenceError(const ConvergenceError &);
    virtual ~ConvergenceError();

private:

    // Not implemented.
    ConvergenceError();
};

#endif // CLASSIC_ConvergenceError_HH
