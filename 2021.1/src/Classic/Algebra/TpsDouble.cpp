// ------------------------------------------------------------------------
// $RCSfile: TpsDouble.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definitions:
//   typedef Tps<double> TpsDouble
//   typedef Vps<double> VpsDouble
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


// Force instantiation of Tps<double> class.
// ------------------------------------------------------------------------

/// Truncated power series in n double variables.
template class Tps<double>;


// Force instantiation of Vps<double> class.
// ------------------------------------------------------------------------

/// Vector of truncated Taylor series in n double variables.
template class Vps<double>;
