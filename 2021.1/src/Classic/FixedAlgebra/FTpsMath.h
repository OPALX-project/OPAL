#ifndef CLASSIC_FTpsMath_HH
#define CLASSIC_FTpsMath_HH

// ------------------------------------------------------------------------
// $RCSfile: FTpsMath.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.4 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// Declared template functions:
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:04:21 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FTps.h"
#include "Utilities/DomainError.h"

#include <complex> // std::real

#include <type_traits>

// 25. March 2017,
// http://stackoverflow.com/questions/30736951/templated-class-check-if-complex-at-compile-time
template<class T> struct is_complex : std::false_type {};
template<class T> struct is_complex<std::complex<T>> : std::true_type {};


// Class FTps; global functions acting on FTps objects.
//   Elementary functions acting on FTps<T,N> objects.
// ------------------------------------------------------------------------

/// Tps x to the power (int y).
template <class T, int N>
FTps<T, N> pow(const FTps<T, N> &x, int y, int trunc = (FTps<T, N>::EXACT));

/// Square root.
template <class T, int N>
FTps<T, N> sqrt(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Sine.
template <class T, int N>
FTps<T, N> sin(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Cosine.
template <class T, int N>
FTps<T, N> cos(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Tangent.
template <class T, int N>
FTps<T, N> tan(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Cotangent.
template <class T, int N>
FTps<T, N> cot(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Secant.
template <class T, int N>
FTps<T, N> sec(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Cosecant.
template <class T, int N>
FTps<T, N> csc(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Exponential.
template <class T, int N>
FTps<T, N> exp(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Natural logarithm.
template <class T, int N>
FTps<T, N> log(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Hyperbolic sine.
template <class T, int N>
FTps<T, N> sinh(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Hyperbolic cosine.
template <class T, int N>
FTps<T, N> cosh(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Hyperbolic tangent.
template <class T, int N>
FTps<T, N> tanh(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Hyperbolic cotangent.
template <class T, int N>
FTps<T, N> coth(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Hyperbolic secant.
template <class T, int N>
FTps<T, N> sech(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Hyperbolic cosecant.
template <class T, int N>
FTps<T, N> csch(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Error function.
template <class T, int N>
FTps<T, N> erf(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));

/// Complementary error function.
template <class T, int N>
FTps<T, N> erfc(const FTps<T, N> &x, int trunc = (FTps<T, N>::EXACT));


// Implementation
// ------------------------------------------------------------------------

template <class T, int N>
FTps<T, N> pow(const FTps<T, N> &x, int y, int trunc) {
    // Default: trunc = EXACT

    FTps<T, N> z(T(1));

    if(y > 0) {
        while(y-- > 0) z = z.multiply(x, trunc);
    } else if(y < 0) {
        if(x[0] == T(0) || x.getMinOrder() != 0)
            throw DomainError("pow(const FTps &,int,int)");
        FTps<T, N> t = x.inverse(trunc);
        while(y++ < 0) z = z.multiply(t, trunc);
    }

    return z;
}


template <class T, int N>
FTps<T, N> sqrt(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::sqrt(FTps<T,N> &x, int trunc)",
                           "Square-root of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if( ( std::real(aZero) <= std::real(T(0)) &&
          std::imag(aZero) <= std::real(T(0))
        ) || x.getMinOrder() != 0)
    {
        std::cerr << "FTps::sqrt(x) called with\nconstant term = " << aZero
                  << "; minOrder = " << x.getMinOrder() << std::endl
                  << "x = " << x << std::endl;
        throw DomainError("sqrt(const FTps &,int)");
    }

    T two_aZero = T(2) * aZero;
    Array1D<T> series(trcOrder + 1);
    series[0] = sqrt(aZero);
    for(int i = 1; i <= trcOrder; i++) {
        series[i] = series[i-1] * double(3 - 2 * i) / (two_aZero * double(i));
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> sin(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::sin(FTps<T,N> &x, int trunc)",
                           "Sine of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(x.getMinOrder() != 0) aZero = T(0);

    Array1D<T> series(trcOrder + 1);
    series[0] = sin(aZero);
    series[1] = cos(aZero);
    for(int i = 2; i <= trcOrder; i++) {
        series[i] = - series[i-2] / double(i * (i - 1));
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> cos(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::cos(FTps<T,N> &x, int trunc)",
                           "Cosine of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(x.getMinOrder() != 0) aZero = T(0);

    Array1D<T> series(trcOrder + 1);
    series[0] =   cos(aZero);
    series[1] = - sin(aZero);
    for(int i = 2; i <= trcOrder; i++) {
        series[i] = - series[i-2] / double(i * (i - 1));
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> tan(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    // The direct series expansion requires the Bernoulli numbers to arbitrary order.
    return sin(x, trunc) / cos(x, trunc);
}

template <class T, int N>
FTps<T, N> cot(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    if(x[0] == T(0) || x.getMinOrder() != 0)
        throw DomainError("cot(const FTps &,int)");

    return cos(x, trunc) / sin(x, trunc);
}


template <class T, int N>
FTps<T, N> sec(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    return cos(x, trunc).inverse();
}


template <class T, int N>
FTps<T, N> csc(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    if(x[0] == T(0) || x.getMinOrder() != 0)
        throw DomainError("csc(const FTps &,int)");

    return sin(x, trunc).inverse();
}


template <class T, int N>
FTps<T, N> exp(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::exp(FTps<T,N> &x, int trunc)",
                           "Exponential of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(x.getMinOrder() != 0) aZero = T(0);

    Array1D<T> series(trcOrder + 1);
    series[0] = exp(aZero);
    for(int i = 1; i <= trcOrder; i++) {
        series[i] = series[i-1] / double(i);
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> log(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::log(FTps<T,N> &x, int trunc)",
                           "Logarithm of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(aZero <= T(0) || x.getMinOrder() != 0)
        throw DomainError("log(const FTps &,int)");

    T a0inv = T(1) / aZero;
    T ain = a0inv;
    Array1D<T> series(trcOrder + 1);
    series[0] = log(aZero);
    series[1] = a0inv;
    for(int i = 2; i <= trcOrder; i++) {
        ain *= -a0inv;
        series[i] = ain / double(i);
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> sinh(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::log(FTps<T,N> &x, int trunc)",
                           "Hyperbolic sine of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(x.getMinOrder() != 0) aZero = T(0);

    Array1D<T> series(trcOrder + 1);
    series[0] = sinh(aZero);
    series[1] = cosh(aZero);
    for(int i = 2; i <= trcOrder; i++) {
        series[i] = series[i-2] / double(i * (i - 1));
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> cosh(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::cosh(FTps<T,N> &x, int trunc)",
                           "Hyperbolic cosine of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(x.getMinOrder() != 0) aZero = T(0);

    Array1D<T> series(trcOrder + 1);
    series[0] = cosh(aZero);
    series[1] = sinh(aZero);
    for(int i = 2; i <= trcOrder; i++) {
        series[i] = series[i-2] / double(i * (i - 1));
    }

    return x.taylor(series, trcOrder);
}


template <class T, int N>
FTps<T, N> tanh(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    return sinh(x, trunc) / cosh(x, trunc);
}


template <class T, int N>
FTps<T, N> coth(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    if(x[0] == T(0) || x.getMinOrder() != 0)
        throw DomainError("coth(const FTps &,int)");

    return cosh(x, trunc) / sinh(x, trunc);
}


template <class T, int N>
FTps<T, N> sech(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    return cosh(x, trunc).inverse();
}


template <class T, int N>
FTps<T, N> csch(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT

    if(x[0] == T(0) || x.getMinOrder() != 0)
        throw DomainError("csch(const FTps &,int)");

    return sinh(x, trunc).inverse();
}

/// Error function.
template <class T, int N>
FTps<T, N> erf(const FTps<T, N> &x, int trunc) {
    
    if ( is_complex<T>::value )
        throw LogicalError("::erf(FTps<T,N> &x, int trunc)",
                           "Error function does not support complex numbers.");
    
    // Default: trunc = EXACT
    int trcOrder = std::min(x.getTruncOrder(), trunc);

    if(trcOrder == FTps<T, N>::EXACT)
        throw LogicalError("::erf(FTps<T,N> &x, int trunc)",
                           "Error function of EXACT polynomial must be truncated.");

    T aZero = x[0];
    if(x.getMinOrder() != 0) aZero = T(0);
    
    Array1D<T> series(trcOrder + 1);
    series[0] = std::erf(std::real(aZero));
    series[1] = 2.0 / std::sqrt(M_PI) * std::exp(-aZero*aZero);
    
    for(int i = 2; i <= trcOrder; ++i) {
        series[i] = - 2.0 / double(i-1) * double((i-2)) * series[i-2] / double(i);
    }
    
    return x.taylor(series, trcOrder);
}

/// Complementary error function.
template <class T, int N>
FTps<T, N> erfc(const FTps<T, N> &x, int trunc) {
    // Default: trunc = EXACT
    
    return T(1) - erf(x, trunc);
}

#endif // CLASSIC_FTpsMath_HH
