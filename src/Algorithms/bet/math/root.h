/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 09-03-2006

   root finding routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#ifndef _BET_ROOT_H
#define _BET_ROOT_H

/* findRoot()
   Originally rtsafe()
   Using a combination of Newton-Raphson and bisection, find the root of a
   function bracketed between x1 and x2. The root, returned as the function
   value rtsafe, will be refined until its accuracy is known within ±xacc.
   funcd is a user-supplied routine that returns both the function value
   and the first derivative of the function. */
double findRoot
(
    void (*funcd)(double, double *, double *), // function to find root
    double, double,                            // boundaries
    double);                                   // accuracy

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
