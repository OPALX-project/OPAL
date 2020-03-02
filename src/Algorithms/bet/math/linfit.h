/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 09-03-2006

   linear fitting routine

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#ifndef _BET_LINFIT_H
#define _BET_LINFIT_H


/* linfit() Given a set of data points x[0..ndata-1],y[0..ndata-1] with
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
            double &sigb, double &chi2, double &q);


/* linfit()
   Given a set of data points x[0..ndata-1],y[0..ndata-1], fit them to
   a straight line y = a + bx by minimizing chi2. Returned are a,b and
   their respective probable uncertainties siga and sigb, and the
   chi-square chi2.
*/
void linfit(double x[], double y[], int ndata,
            double &a, double &b, double &siga,
            double &sigb, double &chi2);

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
