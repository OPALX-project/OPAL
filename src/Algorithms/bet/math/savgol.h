/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 23-11-2006

   Savitzky-Golay Smoothing Filters

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#ifndef _BET_SAVGOL_H
#define _BET_SAVGOL_H

/* sgSmooth()
   Smoothes c[0..n-1] using a Savitzky-Golay filter.
   nl is the number of leftward (past) data points used, while
   nr is the number of rightward (future) data points, making
      the total number of data points used nl+nr+1.
   ld is the order of the derivative desired
      (e.g., ld = 0 for smoothed function).
   m  is the order of the smoothing polynomial,
      also equal to the highest conserved moment;
      usual values are m = 2 or m = 4.
*/
void sgSmooth(
    double [],          // c (data array)
    int,               // number of points [0..n-1]
    int, int,          // nl, nr,
    int, int);         // ld, m

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
