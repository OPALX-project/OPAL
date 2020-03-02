/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 09-03-2006

   root finding routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#include <cmath>

#include "Algorithms/bet/math/root.h"
#include "Utilities/OpalException.h"

/*
   ======================================================================= */
#define MAXIT 100

/* findRoot()
   Originally rtsafe()
   Using a combination of Newton-Raphson and bisection, find the root of a
   function bracketed between x1 and x2. The root, returned as the function
   value rtsafe, will be refined until its accuracy is known within ±xacc.
   funcd is a user-supplied routine that returns both the function value
   and the first derivative of the function. */
double findRoot(
    void (*funcd)(double, double *, double *), // function to find root
    double x1, double x2,                      // boundaries
    double xacc                                // accuracy
    ) {

    double df = 0.0, dx = 0.0, fh =0.0, fl = 0.0;
    (*funcd)(x1, &fl, &df);
    (*funcd)(x2, &fh, &df);

    if((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
        throw OpalException(
            "BetMath_root::findRoot()",
            "Root must be bracketed.");
    }

    if(fl == 0.0) return x1;
    if(fh == 0.0) return x2;

    auto xh = x1;
    auto xl = x2;
    if(fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    auto rts = 0.5 * (x1 + x2);
    auto dxold = std::abs(x2 - x1);
    dx = dxold;
    auto  f = 0.0;
    (*funcd)(rts, &f, &df);
    for(auto j = 1; j <= MAXIT; j++) {
        if((((rts - xh)*df - f) * ((rts - xl)*df - f) >= 0.0)
           || (std::abs(2.0 * f) > std::abs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if(xl == rts)
                return rts;
        } else {
            dxold = dx;
            dx = f / df;
            auto temp = rts;
            rts -= dx;
            if(temp == rts)
                return rts;
        }
        if(std::abs(dx) < xacc)
            return rts;
        (*funcd)(rts, &f, &df);
        if(f < 0.0)
            xl = rts;
        else
            xh = rts;
    }
    throw OpalException(
        "BetMath_root::findRoot()",
        "Maximum number of iterations exeeded.");

    return 0.0;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
