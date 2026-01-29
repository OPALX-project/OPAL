#ifndef CLASSIC_ComplexErrorFun_HH
#define CLASSIC_ComplexErrorFun_HH

// ------------------------------------------------------------------------
// $RCSfile: ComplexErrorFun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Function: Werrf(complex<double>)
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <complex>


// ------------------------------------------------------------------------
/// Complex error function.
//  The algorithms is based on:
//
//  Walter Gautschi, Efficient Computation of the Complex Error Function,
//  [CENTER]
//    SIAM J. Numer. Anal., Vol 7, No. 1, March 1970, pp. 187-198.
//  [/CENTER]
//  Argument: the complex argument,
//  Return value: the complex value of the error function.

extern std::complex<double> Werrf(std::complex<double> z);

#endif // CLASSIC_ComplexErrorFun_HH
