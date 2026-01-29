#ifndef __TUNE__
#define __TUNE__
/*****************************************************************************/
/*                                                                           */
/* TUNE's Class Header                                                       */
/* =====================                                                     */
/*                                                                           */
/*****************************************************************************/

class TUNE_class
/*---------------------------------------------------------------------------*/
{
private:

    // lomb stuff:
    double ofac, hifac;
    double Qmin, Qmax;

public:

    TUNE_class(); //constructor
    virtual ~TUNE_class(void);       //destructor

    int lombAnalysis(double *x, double *y, int Ndat, int nhis);
    int lombAnalysis(std::vector<double> &x, std::vector<double> &y, int nhis, double Norm);


};
/*---------------------------------------------------------------------------*/

#endif
