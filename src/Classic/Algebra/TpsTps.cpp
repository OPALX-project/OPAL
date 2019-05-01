// ------------------------------------------------------------------------
// $RCSfile: TpsTps.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definitions:
//   typedef Tps< Tps<double> > TpsTps
//   typedef Vps< Tps<double> > VpsTps
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


// Force instantiation of Tps< Tps<T> > class.
// ------------------------------------------------------------------------

/// Truncated power series in n Tps<double> variables.
template class Tps< Tps<double> >;


// Force instantiation of Vps<double> class.
// ------------------------------------------------------------------------

/// Vector of truncated power series in n Tps<double> variables.
template class Vps< Tps<double> >;
