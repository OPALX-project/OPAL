/*****************************************************************************/
/*                                                                           */
/* Class for Lomb Periodograms & Co.                                         */
/* =================================                                         */
/*                                                                           */
/* (update: ASM, September 2001)                                             */
/*****************************************************************************/
#include "lomb.h"

#include "Physics/Physics.h"

#include <iostream>

LOMB_class::LOMB_class(int)
/*---------------------------------------------------------------------------*
 * constructor
 * ===========
 *
 *---------------------------------------------------------------------------*/
{}


LOMB_class::~LOMB_class(void)
/*---------------------------------------------------------------------------*
 * destructor
 * ==========
 *
 *---------------------------------------------------------------------------*/
{
    // Do nothing......
}


int LOMB_class::period(std::vector<LOMB_TYPE> *indata, std::vector<LOMB_TYPE> *outdata,
                       double ofac, double hifac, int *nout, int *jmax, double *prob,
                       int amp)
/*---------------------------------------------------------------------------*
 * NR routine
 * ==========
 *
 *---------------------------------------------------------------------------*/
{
    int    i, j, ntmp;
    int    n = 0;
    double  ave, c, cc, cwtau, effm, expy, pnow, pymax, s, ss, sumc, sumcy, sums, sumsh, sumsy;
    double  swtau, var, wtau, xave, xdiff, xmax, xmin, yy;
    //double arg,wtemp,*wi,*wpi,*wpr,*wr;
    double arg, wtemp;

    std::vector<double> wi, wpi, wpr, wr;

    LOMB_TYPE pt;

    CI_lt p, q;
    /*---------------------------------------------------------------------------*/

    wi.erase(wi.begin(), wi.end());
    wpi.erase(wpi.begin(), wpi.end());
    wpr.erase(wpr.begin(), wpr.end());
    wr.erase(wr.begin(), wr.end());

    p = indata->begin();
    q = indata->end();

    n = q - p;

    ntmp = (int)(0.5 * ofac * hifac * n);

    *nout = ntmp;


    if(avevar(indata, &ave, &var) != 0) {
        std::cerr << "LOMB: Average failed!\n";
        return(-1);
    }

    p = indata->begin();
    xmax = xmin = (*p).x;

    for(p = indata->begin(); p != indata->end(); p++) {
        if((*p).x > xmax) xmax = (*p).x;
        if((*p).x < xmin) xmin = (*p).x;
    }

    xdiff = xmax - xmin;
    xave  = 0.5 * (xmax + xmin);
    pymax = 0.0;
    pnow  = 1. / (xdiff * ofac);

    for(p = indata->begin(); p != indata->end(); p++) {

        arg    = Physics::two_pi * (((*p).x - xave) * pnow);
        wpr.push_back((-2. * pow(sin(0.5 * arg) , 2)));
        wpi.push_back(sin(arg));
        wr.push_back(cos(arg));
        wi.push_back(sin(arg));

    }

    // check wr range and data range !!!!
    if((wr.end() - wr.begin()) != n) {
        std::cerr << "LOMB: Vector range mismatch!!!\n";
        return(-1);
    }

    for(i = 0; i < (*nout); i++) {

        pt.x = pnow;
        sumsh = 0. ;
        sumc  = 0.;

        for(j = 0; j < n; j++) {
            c = wr[j];
            s = wi[j];

            sumsh += s * c;
            sumc  += (c - s) * (c + s);
        }

        wtau = 0.5 * atan2(2.0 * sumsh, sumc);
        swtau = sin(wtau);
        cwtau = cos(wtau);
        sums = 0.;
        sumc = 0.;
        sumsy = 0.;
        sumcy = 0.;

        for(j = 0, p = indata->begin() ; j < n && p != indata->end() ; j++, p++) {

            s = wi[j];
            c = wr[j];
            ss = s * cwtau - c * swtau;
            cc = c * cwtau + s * swtau;
            sums += ss * ss;
            sumc += cc * cc;
            yy = (*p).y - ave;
            sumsy += yy * ss;
            sumcy += yy * cc;
            wtemp = wr[j];
            wr[j] = (wr[j] * wpr[j] - wi[j] * wpi[j]) + wr[j];
            wi[j] = (wi[j] * wpr[j] + wtemp * wpi[j]) + wi[j];

        }

        pt.y = 0.5 * (sumcy * sumcy / sumc + sumsy * sumsy / sums) / var;
        if(amp) pt.y = sqrt(pow((sumcy / sumc), 2.) + pow((sumsy / sums), 2.));

        if(pt.y >= pymax) {
            pymax = pt.y;
            *jmax = i;
        }

        outdata->push_back(pt);

        pnow += 1. / (ofac * xdiff);

    }

    expy  = exp(-pymax);
    effm  = 2. * (*nout) / ofac;
    *prob = effm * expy;

    if(*prob > 0.01) *prob = 1. - pow((1. - expy), effm);

    wi.erase(wi.begin(), wi.end());
    wpi.erase(wpi.begin(), wpi.end());
    wpr.erase(wpr.begin(), wpr.end());
    wr.erase(wr.begin(), wr.end());

    return(0);

}


int LOMB_class::avevar(std::vector<LOMB_TYPE> *data, double *ave, double *var)
/*---------------------------------------------------------------------------*
 * NR routine
 * ==========
 *
 *---------------------------------------------------------------------------*/
{
    int    n;
    double  s, ep;

    CI_lt p, q;

    /*---------------------------------------------------------------------------*/


    *ave = 0.;
    p = data->begin();
    q = data->end();

    n = q - p;

    if(n < 2) {
        std::cerr << "Only one datapoint -> no averaging....\n";
        return(-1);
    }

    for(p = data->begin(); p != data->end(); p++) *ave += (*p).y;

    *ave = *ave / n;
    *var = 0.;
    ep   = 0.;

    for(p = data->begin(); p != data->end(); p++) {
        s     = (*p).y - *ave;
        ep   += s;
        *var += s * s;
    }

    *var = (*var - ep * ep / n) / (n - 1);

    return(0);

}



double LOMB_class::signi(double *peak, int *nout, double *ofac)
/*---------------------------------------------------------------------------*
 * Calculate the significance of a peak in an Lomb Periodogram
 * ===========================================================
 *
 * Input:  double* peak: Peak of periodogram
 *         int*    nout: Number of frequencies
 *         double* ofac: Oversampling factor
 * Output: double  Sign: Significance of peak
 *---------------------------------------------------------------------------*/
{

    double  expy, effm, prob;

    /*---------------------------------------------------------------------------*/

    expy = exp(-1 * (*peak));
    effm = 2. * (double)(*nout) / (*ofac);
    prob = effm * expy ;
    if(prob > 0.01) prob = 1. - pow((1. - expy), effm);

    return(prob);


}


int LOMB_class::moment(std::vector<LOMB_TYPE> *indata, double *ave, double *adev,
                       double *sdev, double *var, double *skew, double *curt)
/*---------------------------------------------------------------------------*
 * Calculate the first moments of a distribution (free after NR)
 * =============================================================
 *
 * Input:  vector indata: (periodogram) value vector of type LOMB_TYPE
 * Output: double* ave  : average
 *         double* adev : average deviation
 *         double* sdev : standard deviation
 *         double* var  : variance
 *         double* skew : skewness
 *         double* curt : kurtosis
 *---------------------------------------------------------------------------*/
{

    int      n;
    double  pnr, s, ep;

    CI_lt p, q;
    /*---------------------------------------------------------------------------*/

    p = indata->begin();
    q = indata->end();

    n = q - p;

    if(n < 2) {
        std::cerr << "To few data points for moment analysis!\n";
        return(-1);
    }


    /*
     * First pass to get the mean
     * --------------------------
     */

    s = 0;
    double nn = 0;
    for(p = indata->begin(); p != indata->end(); p++) {
        s += (*p).y * (*p).x;
        nn += (*p).y;
    }

    *ave = s / (double)nn;

    /*
     * Second  pass
     * ------------
     */


    *adev = 0.;
    *var  = 0.;
    *skew = 0.;
    *curt = 0.;
    ep    = 0.;

    for(p = indata->begin(); p != indata->end(); p++) {

        s = ((*p).x - *ave);
        ep += s * (*p).y;
        *adev = *adev + (double)std::abs((double)s);

        pnr = s * s;
        *var += pnr * (*p).y;

        pnr *= s;
        *skew += pnr * (*p).y;

        pnr *= s;
        *curt += pnr * (*p).y;

    }

    *adev = *adev / (double)nn;

    *var = (*var - ep * ep / (double)nn) / ((double)(nn - 1));

    *sdev = sqrt(*var);

    if(*var != 0.) {
        *skew = *skew / ((double)nn * pow(*sdev, 3.));
        *curt = *curt / ((double)nn * pow(*var, 2.)) - 3.;
    } else {
        std::cerr << "No skew or kurtosis when zero variance in moment\n";
    }

    return(0);

}
