/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 07-03-2006

   Runge-Kutta with adaptive stepsize control
   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#ifndef _BET_RK_H
#define _BET_RK_H

/*
   Given values for the variables y[0..n-1] use the fourth-order Runge-Kutta
   method to advance the solution over an interval h and return the
   incremented variables in y[0..n-1]. The user supplies the routine
   derivs(x,y,dydx), which returns derivatives dydx at x.t
*/
void rk4(
    double y[],    // initial condition
    int n,     // number of equations
    double x,      // start
    double h,      // interval
    void (*derivs)(double, double [], double []));

/*
   Runge-Kutta driver with adaptive stepsize control.

   Integrate starting values ystart[0..nvar-1]
   from x1 to x2 with accuracy eps,

   storing intermediate results in global variables.
   h1 should be set as a guessed first stepsize,
   hmin as the minimum allowed stepsize (can be zero).

   On output:
   - nok and nbad are the number of good and bad (but retried and fixed)
     steps taken
   - ystart is replaced by values at the end of the integration interval.
   - RETURNS 0 on success and 1 on failure

   derivs is the user-supplied routine for calculating the right-hand side
   derivative.

   rkqs is the name of the stepper routine (see above).
*/
bool odeint(
    double ystart[],  // initial condition
    int    nvar,      // number of equations
    double x1,        // start
    double x2,        // end
    double eps,       // accuracy
    double h1,        // guessed first step-size
    double hmin,      // minimum step-size
    int    &nok,      // number of good steps
    int    &nbad,     // number of bad steps
    void (*derivs)(double, double [], double []));

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
