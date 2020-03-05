/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 07-03-2006

   Runge-Kutta with adaptive stepsize control

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#include <cmath>

#include "Utility/IpplInfo.h"
#include "Algorithms/bet/math/rk.h"

/* Internal functions
   ========================================================================= */

/*
   Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1]
   known at x, use the fifth-order Cash-Karp Runge-Kutta method
   to advance the solution over an interval h and return the incremented
   variables as yout[0..n-1].

   Also return an estimate of the local truncation error in yout
   using the embedded fourth-order method.

   The user supplies the routine derivs(x,y,dydx),
   which returns derivatives dydx at x.
*/
static void rkck (
    double y[],
    double dydx[],
    int    n,
    double x,
    double h,
    double yout[],
    double yerr[],
    void (*derivs)(double, double [], double [])
    ) {

    const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2;
    const double b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2;
    const double b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0;
    const double b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0;
    const double b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0;
    const double c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0;
    const double dc5 = -277.0 / 14336.0;
    const double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0;
    const double dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;
    double ak2[n];
    double ak3[n];
    double ak4[n];
    double ak5[n];
    double ak6[n];
    double ytemp[n];

    for(auto i = 0; i < n; i++)
        ytemp[i] = y[i] + b21 * h * dydx[i];
    (*derivs)(x + a2 * h, ytemp, ak2);
    for(auto i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    (*derivs)(x + a3 * h, ytemp, ak3);
    for(auto i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    (*derivs)(x + a4 * h, ytemp, ak4);
    for(auto i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    (*derivs)(x + a5 * h, ytemp, ak5);
    for(auto i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    (*derivs)(x + a6 * h, ytemp, ak6);
    for(auto i = 0; i < n; i++)
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
    for(auto i = 0; i < n; i++)
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
}


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
// The value ERRCON equals (5/SAFETY) raised to the power (1/PGROW), see use below.

/*
   Fifth-order Runge-Kutta step with monitoring of local truncation error
   to ensure accuracy and adjust stepsize.
   Input are:
   - the dependent variable vector y[0..n-1]
   - its derivative dydx[0..n-1]
   at the starting value of the independent variable x.

   Also input are:
   - the stepsize to be attempted htry
   - the required accuracy eps
   - the vector yscal[0..n-1] against which the error is scaled.

   On output, y and x are replaced by their new values,
   hdid is the stepsize that was actually accomplished,
   and hnext is the estimated next stepsize.

   derivs is the user-supplied routine that computes the
   right-hand side derivatives.
*/
static bool rkqs
(
    double y[], double dydx[], int n, double &x, double htry, double eps,
    double yscal[], double &hdid, double &hnext,
    void (*derivs)(double, double [], double [])
    ) {

    double yerr[n];
    double ytemp[n];;
    double h = htry;
    int nLoop = 0;
    for(;;) {
        ++nLoop;
        rkck(y, dydx, n, x, h, ytemp, yerr, derivs);
        double errmax = 0.0;
        for(auto i = 0; i < n; i++)
            errmax = std::max(errmax, std::abs(yerr[i] / yscal[i]));
        errmax /= eps;
        if(errmax > 1.0) {
            h = SAFETY * h * std::pow(errmax, PSHRNK);
            if(h < 0.1 * h)
                h *= 0.1;
            auto xnew = x + h;
            if(xnew == x) {
                WARNMSG("BetMath_rk::rkqs(): stepsize underflow");
                return false;
            }
            continue;
        } else {
            if(errmax > ERRCON)
                hnext = SAFETY * h * std::pow(errmax, PGROW);
            else
                hnext = 5.0 * h;
            x += (hdid = h);
            for(auto i = 0; i < n; i++)
                y[i] = ytemp[i];
            break;
        }
    }
    return true;
}


/* External functions
   ========================================================================= */

/*
   Given values for the variables y[0..n-1] use the fourth-order Runge-Kutta
   method to advance the solution over an interval h and return the
   incremented variables in y[0..n-1]. The user supplies the routine
   derivs(x,y,dydx), which returns derivatives dydx at x.t
*/
void rk4(double y[], int n, double x, double h,
         void (*derivs)(double, double [], double [])) {

    double dydx[n];
    double dym[n];
    double dyt[n];
    double yt[n];

    (*derivs)(x, y, dydx);
    auto hh = h * 0.5;
    auto h6 = h / 6.0;
    auto xh = x + hh;
    for(auto i = 0; i < n; i++)
        yt[i] = y[i] + hh * dydx[i];
    (*derivs)(xh, yt, dyt);
    for(auto i = 0; i < n; i++)
        yt[i] = y[i] + hh * dyt[i];
    (*derivs)(xh, yt, dym);
    for(auto i = 0; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
        dym[i] += dyt[i];
    }
    (*derivs)(x + h, yt, dyt);
    for(auto i = 0; i < n; i++)
        y[i] += h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
}

/*
   Runge-Kutta driver with adaptive stepsize control.

   Integrate starting values ystart[0..nvar-1]
   from x1 to x2 with accuracy eps,

   On output:
   - nok and nbad are the number of good and bad (but retried and fixed)
     steps taken
   - ystart is replaced by values at the end of the integration interval.
   - RETURNS 0 on success and 1 on error

   derivs is the user-supplied routine for calculating the right-hand side
   derivative.

   rkqs is the name of the stepper routine (see above).
*/
#define MAXSTP 10000
#define TINY 1.0e-30

bool odeint(
    double ystart[],  // initial condition
    int    nvar,      // number of equations
    double x1,
    double x2,
    double eps,
    double h1,
    double hmin,
    int    &nok,
    int    &nbad,
    void (*derivs)(double, double [], double [])
    ) {

    double yscal[nvar];
    double y[nvar];
    double dydx[nvar];
    auto x = x1;
    auto h = ((x2-x1) >= 0.0 ? std::abs(h1) : -std::abs(h1));
    nok = nbad = 0;
    for(auto i = 0; i < nvar; i++)
        y[i] = ystart[i];
    for(auto nstp = 0; nstp < MAXSTP; nstp++) {
        (*derivs)(x, y, dydx);
        for(auto i = 0; i < nvar; i++) {
            yscal[i] = std::abs(y[i]) + std::abs(dydx[i] * h) + TINY;
        }
        if((x + h - x2) * (x + h - x1) > 0.0)
            h = x2 - x;
        double hdid = 0.0, hnext = 0.0;
        rkqs(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs);
        if(hdid == h)
            ++nok;
        else
            ++nbad;
        if((x - x2) * (x2 - x1) >= 0.0) {
            for(auto i = 0; i < nvar; i++)
                ystart[i] = y[i];
            return true;
        }
        if(std::abs(hnext) <= hmin) {
            WARNMSG("BetMath_rk::odeint(): Step size too small.");
            return false;
        }
        h = hnext;
    }
    WARNMSG("BetMath_rk::odeint(): Too many steps in routine.");
    return false;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
