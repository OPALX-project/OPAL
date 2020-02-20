/*****************************************************************************/
/*                                                                           */
/* Class TUNE                                                                */
/* ==============                                                            */
/*                                                                           */
/* ASM, September 2001                                                       */
/*****************************************************************************/
#include <algorithm>
#include <vector>
#include <cstring>

#include "Utility/Inform.h"

#include "Algorithms/Ctunes.h"
#include "Algorithms/lomb.h"

extern Inform *gmsg;

//RANLIB_class rndm(265314159,4);


TUNE_class::TUNE_class():
        ofac(0.0),
        hifac(0.0),
        Qmin(0.0),
        Qmax(0.0)
/*---------------------------------------------------------------------------*
 * constructor
 * ===========
 *
 *---------------------------------------------------------------------------*/
{
}

TUNE_class::~TUNE_class(void)
/*---------------------------------------------------------------------------*
 * destructor
 * ==========
 *
 *---------------------------------------------------------------------------*/
{

    // Do nothing......

}

int TUNE_class::lombAnalysis(std::vector<double> &x, std::vector<double> &y, int /*nhis*/, double Norm)
/*-----------------------------------------------------------------------------
 *  Launch Lomb analysis and plot results
 *  =======================================
 *
 *---------------------------------------------------------------------------*/
{

    int Ndat = x.size();
    int    i, nout, jmax;
    int    pairc;
    int    datcnt = 0;
    int stat = 0;
    double prob, probi;
    double tofac = 0.8;

    LOMB_TYPE tlom;

    CI_lt p, q, r, s, tp;

    char   mess[80];
    LOMB_class *la;
    std::vector<LOMB_TYPE> lodata, lodata2;
    /*---------------------------------------------------------------------------*/

    /*
     * Do Lomb analysis
     * ================
     */

    for(int j = 0; j < Ndat; j++) {
        tlom.x = x[j];
        tlom.y = y[j];
        lodata.push_back(tlom);
    }

    p = lodata.begin();
    q = lodata.end();

    datcnt = (int) count_if(p, q, Lomb_eq(0.));

    if(datcnt > (q - p - 10)) {
        memset(mess, '\0', sizeof(mess));
        sprintf(mess, "Just found %d data points that are == 0!", datcnt);
        *gmsg << "* " << mess << endl;
        return(-1);
    }

    // this parameterset works ok in most cases.....
    ofac  = 4.0;
    hifac = 0.8;
    Qmin  = 0.2;
    Qmax  = 0.4;

    la = new LOMB_class(1);

    stat = 0;
    stat = la->period(&lodata, &lodata2, ofac, hifac, &nout, &jmax, &prob, 0);
    if(stat != 0) {
        memset(mess, '\0', sizeof(mess));
        sprintf(mess, "@C3ERROR: Lomb analysis failed!");
        *gmsg << "* " << mess << endl;

        delete la;
        la = NULL;
        return(-1);
    }

    double pairx[nout];
    double pairy[nout];

    pairc = 0;
    for(i = 0; i < nout; i++) {
        if(lodata2[i].y > 2.) {
            pairx[pairc] = lodata2[i].x;
            pairy[pairc] = lodata2[i].y;
            if((pairy[pairc] > pairy[pairc-1]) &&
               (pairy[pairc] > lodata2[i+1].y)) {
                probi = la->signi(&pairy[pairc], &nout, &tofac);
                if(pairy[pairc] > 4.) {
                    memset(mess, '\0', sizeof(mess));
                    sprintf(mess, "%12.8f %8.2f %8.3f %d", pairx[pairc]*Norm, pairy[pairc], probi, i);
                    *gmsg << "* " << mess << endl;
                }
            }
            pairc++;
        }
    }

    memset(mess, '\0', sizeof(mess));
    sprintf(mess, " ===> Max: %12.8f %8.3f\n", lodata2[jmax].x * Norm, lodata2[jmax].y);
    *gmsg << "* " << mess << endl;

    delete la;
    la = NULL;
    return(0);
}


int TUNE_class::lombAnalysis(double *x, double *y, int Ndat, int /*nhis*/)
/*-----------------------------------------------------------------------------
 *  Launch Lomb analysis and plot results
 *  =======================================
 *
 *---------------------------------------------------------------------------*/
{
    int    i, nout, jmax;
    int    pairc;
    int    datcnt = 0;
    int stat = 0;
    double prob, probi;
    double tofac = 0.8;

    LOMB_TYPE tlom;

    CI_lt p, q, r, s, tp;

    char   mess[80];
    LOMB_class *la;
    std::vector<LOMB_TYPE> lodata, lodata2;
    /*---------------------------------------------------------------------------*/

    sprintf(mess, "TUNE_class LombAnalysis requested");
    *gmsg << "* " << mess << endl;

    /*
     * Do Lomb analysis
     * ================
     */

    for(int j = 0; j < Ndat; j++) {
        tlom.x = x[j];
        tlom.y = y[j];
        lodata.push_back(tlom);
    }

    p = lodata.begin();
    q = lodata.end();

    datcnt = count_if(p, q, Lomb_eq(0.));

    if(datcnt > (q - p - 10)) {
        memset(mess, '\0', sizeof(mess));
        sprintf(mess, "Just found %d data points that are == 0!", datcnt);
        *gmsg << "* " << mess << endl;
        return(-1);
    }

    // this parameterset works ok in most cases.....
    ofac  = 4.0;
    hifac = 0.8;
    Qmin  = 0.2;
    Qmax  = 0.4;

    la = new LOMB_class(1);

    stat = 0;
    stat = la->period(&lodata, &lodata2, ofac, hifac, &nout, &jmax, &prob, 0);
    if(stat != 0) {
        memset(mess, '\0', sizeof(mess));
        sprintf(mess, "@C3ERROR: Lomb analysis failed!");
        *gmsg << "* " << mess << endl;

        delete la;
        la = NULL;
        return(-1);
    }

    memset(mess, '\0', sizeof(mess));
    sprintf(mess, "=====> jmax = %d", jmax);
    *gmsg << "* " << mess << endl;

    double pairx[nout];
    double pairy[nout];

    memset(mess, '\0', sizeof(mess));
    sprintf(mess, "\n********** Peaks in Data:       **************");
    *gmsg << "* " << mess << endl;

    /*
    ada make histogram
    hbook1(nhis,"Lomb data",nout,
       (float)lodata2[0].x,
       (float)lodata2[nout-1].x);

    */
    pairc = 0;
    for(i = 0; i < nout; i++) {
        /* ada book histogram
        Hf1(nhis,(float)lodata2[i].x,(float)lodata2[i].y);
        */
        if(lodata2[i].y > 2.) {
            pairx[pairc] = lodata2[i].x;
            pairy[pairc] = lodata2[i].y;
            if((pairy[pairc] > pairy[pairc-1]) &&
               (pairy[pairc] > lodata2[i+1].y)) {
                probi = la->signi(&pairy[pairc], &nout, &tofac);
                if(pairy[pairc] > 4.) {
                    memset(mess, '\0', sizeof(mess));
                    sprintf(mess, "%12.8f %8.2f %8.3f %d", pairx[pairc], pairy[pairc],
                            probi, i);
                    *gmsg << "* " << mess << endl;
                }
            }
            pairc++;
        }
    }

    memset(mess, '\0', sizeof(mess));
    sprintf(mess, "\n===> Max: %12.8f %8.3f\n", lodata2[jmax].x, lodata2[jmax].y);
    *gmsg << "* " << mess << endl;

    delete la;
    la = NULL;

    return(0);

}