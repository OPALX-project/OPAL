/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 09-03-2006

   linear fitting routine

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#include <cmath>

#include "Utilities/OpalException.h"
#include "Algorithms/bet/math/linfit.h"

/* Internal functions
   ========================================================================= */

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

/*
   Returns the ln-value of the gamma function ln[G(xx)] for xx > 0.
*/
static double gammln(double xx) {
    const double cof[6] = {76.18009172947146, -86.50532032941677,
                           24.01409824083091, -1.231739572450155,
                           0.1208650973866179e-2, -0.5395239384953e-5};
    auto x = xx;
    auto y = xx;
    auto tmp  = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser  = 1.000000000190015;
    for(auto j = 0; j <= 5; j++)
        ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}


/*
   Returns the incomplete gamma function Q(a, x) evaluated by its
   continued fraction representation as gammcf. Also returns Gamma(a) as
   gln.
*/

static double gcf(double a, double x) {

    auto gln = gammln(a);
    auto b = x + 1.0 - a;
    auto c = 1.0 / FPMIN;
    auto d = 1.0 / b;
    auto h = d;
    int i = 1;
    while(i <= ITMAX) {
        auto an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if(std::abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if(std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        auto del = d * c;
        h *= del;
        if(std::abs(del - 1.0) < EPS) break;
        i++;
    }
    if(i > ITMAX) {
        throw OpalException(
            "BetMath_linfit::gcf()",
            "linfit-->gfc");
    }
    return exp(-x + a * log(x) - gln) * h;
}


/* 
   Returns the incomplete gamma function P(a,x) evaluated by
   its series representation as gamser.
*/
static double gser(double a, double x) {

    auto gln = gammln(a);
    if(x <= 0.0) {
        if(x < 0.0) {
            throw OpalException (
                "BetMath_linfit::gsr()",
                "linfit-->gser");
        }
        return 0.0;
    } else {
        auto ap = a;
        auto sum = 1.0 / a;
        auto del = 1.0 / a;
        for(auto n = 1; n <= ITMAX; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if(std::abs(del) < std::abs(sum)*EPS) {
                return sum * exp(-x + a * log(x) - gln);
            }
        }
        throw OpalException (
            "BetMath_linfit::gser()",
            "linfit-->gser");
    }
}

/*
   Returns the incomplete gamma function Q(a, x) = 1 - P(a, x).
*/
static double gammq(double a, double x) {

    if(x < 0.0 || a <= 0.0) {
        throw OpalException (
            "BetMath_linfit::gammq()",
            "linfit-->gammg");
    }
    if(x < (a + 1.0)) {
        return 1.0 - gser(a, x);
    } else {
        return gcf(a, x);
    }
}

/* external functions
   ========================================================================= */

/* 
   Given a set of data points x[0..ndata-1],y[0..ndata-1] with
   individual standard deviations sig[0..ndata-1], fit them to a
   straight line y = a + bx by minimizing chi2. Returned are a,b and
   their respective probable uncertainties siga and sigb, the
   chi-square chi2, and the goodness-of-fit probability q (that the
   fit would have chi2 this large or larger). If mwt=0 on input, then
   the standard deviations are assumed to be unavailable: q is
   returned as 1.0 and the normalization of chi2 is to unit standard
   deviation on all points.
*/
void linfit(double x[], double y[], int ndata,
            double sig[], int mwt, double &a, double &b, double &siga,
            double &sigb, double &chi2, double &q) {

    double sx = 0.0, sy = 0.0, st2 = 0.0;

    b = 0.0;
    double ss = 0.0;
    if(mwt) {
        for(auto i = 0; i < ndata; i++) {
            auto wt = 1.0 / std::pow(sig[i], 2);
            ss += wt;
            sx += x[i] * wt;
            sy += y[i] * wt;
        }
    } else {
        for(auto i = 0; i < ndata; i++) {
            sx += x[i];
            sy += y[i];
        }
        ss = ndata;
    }
    auto sxoss = sx / ss;
    if(mwt) {
        for(auto i = 0; i < ndata; i++) {
            auto t = (x[i] - sxoss) / sig[i];
            st2 += t * t;
            b += t * y[i] / sig[i];
        }
    } else {
        for(auto i = 0; i < ndata; i++) {
            auto t = x[i] - sxoss;
            st2 += t * t;
            b += t * y[i];
        }
    }
    b /= st2;
    a = (sy - sx * b) / ss;
    siga = sqrt((1.0 + sx * sx / (ss * st2)) / ss);
    sigb = sqrt(1.0 / st2);
    chi2 = 0.0;
    if(mwt == 0) {
        for(auto i = 0; i < ndata; i++)
            chi2 += std::pow(y[i] - a - b * x[i], 2);
        q = 1.0;
        auto sigdat = sqrt(chi2 / (ndata - 2));
        siga *= sigdat;
        sigb *= sigdat;
    } else {
        for(auto i = 0; i < ndata; i++)
            chi2 += std::pow((y[i] - a - b * x[i]) / sig[i], 2);
        q = gammq(0.5 * (ndata - 2), 0.5 * chi2);
    }
}

/*
   Given a set of data points x[0..ndata-1],y[0..ndata-1] with, fit
   them to a straight line y = a + bx by minimizing chi2. Returned are
   a,b and their respective probable uncertainties siga and sigb, and
   the chi-square chi2.
*/
void linfit(double x[], double y[], int ndata,
            double &a, double &b, double &siga,
            double &sigb, double &chi2) {

    b = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    for(auto i = 0; i < ndata; i++) {
        sx += x[i];
        sy += y[i];
    }

    auto sxoss = sx / ndata;

    double st2 = 0.0;
    for(auto i = 0; i < ndata; i++) {
        auto t = x[i] - sxoss;
        st2 += t * t;
        b += t * y[i];
    }
    b /= st2;
    a = (sy - sx * b) / ndata;
    siga = sqrt((1.0 + sx * sx / (ndata * st2)) / ndata);
    sigb = sqrt(1.0 / st2);
    chi2 = 0.0;

    for(auto i = 0; i < ndata; i++)
        chi2 += std::pow(y[i] - a - b * x[i], 2);

    auto sigdat = std::sqrt(chi2 / (ndata - 2));
    siga *= sigdat;
    sigb *= sigdat;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
