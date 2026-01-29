// ------------------------------------------------------------------------
// $RCSfile: TpsComplex.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definitions:
//   Tps<std::complex<double> > TpsComplex
//   Vps<std::complex<double> > VpsComplex
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Tps.hpp"
#include "Algebra/Vps.hpp"
#include <complex>


// Force instantiation of Tps<std::complex<double> > class.
// ------------------------------------------------------------------------

/// Truncated power series in n complex variables
template class Tps<std::complex<double> >;


// Force instantiation of Vps<std::complex<double> > class.
// ------------------------------------------------------------------------

/// Vector of truncated Taylor series in n complex variables.
template class Vps<std::complex<double> >;
