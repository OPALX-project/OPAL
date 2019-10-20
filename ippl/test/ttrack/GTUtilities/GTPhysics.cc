// ------------------------------------------------------------------------
// $RCSfile: GTPhysics.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Physics:
//   Contains the definitions for the mathematical and physical constants.
//
// ------------------------------------------------------------------------
// Class category: Physics
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "GTUtilities/GTPhysics.h"

// Namespace Physics
// ------------------------------------------------------------------------

namespace Physics {
  
  // Universal mathematical constants
  // ----------------------------------------------------------------------
  const double pi         = 3.14159265358979323846;
  const double two_pi     = 2.0 * pi;
  const double u_two_pi   = 1.0 / (2.0 * pi);
  const double e          = 2.7182818284590452354;
  const double log10e     = 0.43429448190325182765;
  
  
  // Universal physical constants
  // ----------------------------------------------------------------------
  
  // Global constants:
  const double c          = 299792458.0;
  const double mu_0       = 1.25663706212e-06;
  const double epsilon_0  = 8.8541878128e-12;
  const double h_bar      = 6.582119569e-25;
  
  // electromagnetic constants:
  const double q_e        = 1.602176634e-19;
  const double alpha      = 7.2973525693e-03;
  
  // constants for electrons:
  const double m_e        = 0.51099895000e-03;
  const double r_e        = 2.8179403262e-15;
  const double lamda_e    = 3.8615926796e-13;
  const double a_e        = 1.15965218128e-03;
  
  // constants for protons:
  const double m_p        = 0.93827208816e+00;
  const double r_p        = 1.53469857e-18;
  const double lamda_p    = 2.10308910336e-16;
  const double a_p        = 1.792847386e+00;

  const double PMASS       = 1.67262192369e-27;  // kg
  const double EMASS       = 9.1093837015e-31; // kg
  const double PCHARGE     = 1.602176634e-19; // C
};
