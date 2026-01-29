// ------------------------------------------------------------------------
// $RCSfile: ComplexErrorFun.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Function: Werrf(complex<double>)
//   The complex error function.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:33:44 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Utilities/ComplexErrorFun.h"
#include <complex>


// Complex error function
// The algorithms is based on:
//   Walter Gautschi,
//   Efficient Computation of the Complex Error Function,
//   SIAM J. Numer. Anal., Vol 7, No. 1, March 1970, pp. 187-198.
// ------------------------------------------------------------------------

std::complex<double> Werrf(std::complex<double> z) {
    const double cc   = 1.12837916709551;  // 2 / sqrt(pi)
    const double xlim = 5.33;
    const double ylim = 4.29;

    const double x = std::abs(std::real(z));
    const double y = std::abs(std::imag(z));
    std::complex<double> w;

    if(y < ylim  &&  x < xlim) {
        // Inside limit rectangle: equation (3.8): q = s(z);
        const double q = (1.0 - y / ylim) * sqrt(1.0 - (x * x) / (xlim * xlim));
        // Equations (3.11).
        const int nc = 6 + int(23.0 * q);
        const int nu = 9 + int(21.0 * q);
        const std::complex<double> zz(1.6 * q + y, - x); // h - i*z

        // Equations (3.12) for n = nu, nu - 1, ..., N.
        std::complex<double> r(0.0);               // r_{nu}
        for(int n = nu; n > nc; n--) {             // for n = nu, nu-1, ..., nc+1
            r = 0.5 / (zz + double(n + 1) * r);      // r_{n-1}
        }

        // Equations (3.12) for n = N, N - 1, ..., 0.
        double power = pow(3.2 * q, nc);           // (2*h)^nc
        const double div   = 0.3125 / q;           // 1 / (2*h)
        std::complex<double> s(0.0);               // s_{N}
        for(int n = nc; n >= 0; n--) {             // for n = N, N-1, ..., 0
            r = 0.5 / (zz + double(n + 1) * r);      // r_{n-1}
            s = r * (power + s);                     // s_{n-1}
            power *= div;                            // (2*h)^(n-1)
        }

        // Equation (3.13).
        w = cc * s;
    } else {
        // Outside limit rectangle: equations (3.12) for h = 0.
        const std::complex<double> zz(y, - x);     // - i*z
        std::complex<double> r(0.0);               // r_{nu}
        for(int n = 8; n >= 0; n--) {              // for n = nu, nu-1, ..., 0
            r = 0.5 / (zz + double(n + 1) * r);      // r_{n-1}
        }

        // Equation (3.13).
        w = cc * r;
    }

    // Equations (3.1).
    if(std::imag(z) < 0.0) {
        std::complex<double> zz(x, y);
        w = 2.0 * std::exp(- zz * zz) - w;
        if(std::real(z) > 0.0) w = std::conj(w);
    } else {
        if(std::real(z) < 0.0) w = std::conj(w);
    }

    return w;
}
