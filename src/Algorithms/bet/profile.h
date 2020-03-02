/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 07-03-2006

   calculates a functional profile from a mapping
*/

#ifndef _BET_PROFILE_H
#define _BET_PROFILE_H

#include <vector>

#include <gsl/gsl_spline.h>

enum Interpol_type {
    itype_spline,      // spline interpolation
    itype_lin,         // linearl interpolation
};

class Profile {
    int
    n;
    double
    yMax,                // maximum of array
    yMin,                // minimum of array
    sf;                  // scaling factor (1.0 by default)
    std::vector<double> x, y, y2;
    gsl_interp_accel* acc;
    gsl_spline* spline;

private:
    void create();         // general creator routine

public:
    Profile(double);       // dummy creator (returns given value on all requests)
    Profile(               // creator from array
        double *,              // x
        double *,              // y
        int);                  // number of points
    Profile(               // creator from file
        char *,                 // filename
        double = 0.0);          // cutoff value
    ~Profile();
    
    void normalize();      // set max of profile to 1.0
    void scale(double);    // scale the amplitude
    double set(double);    /* set the amplitude
                returns the new scaling factor sf */
    void   setSF(double);  // set sf
    double getSF();        // get sf

    double get(            // get a value
        double,                 // x
        Interpol_type = itype_spline);

    int    getN();               // get number of points

    double max();                // get maximum y-value
    double min();                // get minimum y-value
    double xMax();               // get maximum x-value
    double xMin();               // get minimum x-value

    double Leff();               // get effective length
    double Leff2();              // get effective length of y^2
    double Labs();               // get effective length from absolute value of field

    void   dump(FILE * = stdout, double = 0.0);
    void   dump(          // dump the profile in a SDDS file (ascii)
        char *,                 // filename
        double = 0.0);          // offset

};

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
