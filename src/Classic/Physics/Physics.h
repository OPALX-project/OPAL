#ifndef CLASSIC_Physics_HH
#define CLASSIC_Physics_HH

// ------------------------------------------------------------------------
// $RCSfile: Physics.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Physics
//
// ------------------------------------------------------------------------
// Class category: Physics
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class Physics
// ------------------------------------------------------------------------
/// A namespace defining various mathematical and physical constants.

namespace Physics {

    /// The value of \f[ \pi \f]
    constexpr double pi         = 3.14159265358979323846;

    /// The value of \f[2 \pi \f]
    constexpr double two_pi     = 2 * pi;

    /// The value of \f[ \frac{1}{2} \pi \f]
    constexpr double u_two_pi   = 1.0 / two_pi;

    /// The value of \f[ e \f]
    constexpr double e          = 2.7182818284590452354;

    /// The logarithm of $e$ to the base 10
    constexpr double log10e     = 0.43429448190325182765;

    /// The conversion factor from radians to degrees
    constexpr double rad2deg    = 180.0 / pi;

    /// The conversion factor from degrees to radians
    constexpr double deg2rad    = 1.0 / rad2deg;

    /// The velocity of light in m/s
    constexpr double c          = 299792458.0;

    /// The permeability of vacuum in Vs/Am
    constexpr double mu_0       = 1.256637061e-06;

    /// The permittivity of vacuum in As/Vm
    constexpr double epsilon_0  = 8.854187817e-12;

    /// The reduced Planck constant in GeVs
    constexpr double h_bar      = 6.5821220e-25;

    /// The Avogadro's number
    constexpr double Avo        = 6.022e23;

    /// Boltzman's constant in eV/K.
    constexpr double kB         = 8.6173324e-5;

    /// The elementary charge in As
    constexpr double q_e        = 1.60217733e-19;

    /// The fine structure constant, no dimension
    constexpr double alpha      = 7.29735308e-03;

    /// The atomic mass constant energy equivalent in GeV
    constexpr double uam        = 0.93149410242;

    /// The electron rest mass in GeV
    constexpr double m_e        = 0.51099892e-03;

    /// Rydberg's energy (Rydberg's constant times hc) in GeV
    constexpr double E_ryd      = 13.605693122994e-9;
 
    /// The classical electron radius in m
    constexpr double r_e        = 2.81794092e-15;

    /// Bohr radius in m
    constexpr double a0         = 5.29177210903e-11;

    /// The reduced Compton wave length for electrons in m
    constexpr double lamda_e    = 3.86159323e-13;

    /// The magnetic momentum anomaly for electrons, no dimension
    constexpr double a_e        = 1.159652193e-03;

    /// The proton rest mass in GeV
    constexpr double m_p        = 0.93827204e+00;

    /// The classical proton radius in m
    constexpr double r_p        = 1.53469857e-18;

    /// The reduced Compton wave length for protons in m
    constexpr double lamda_p    = 2.10308937e-16;

    /// The magnetic momentum anomaly for protons, no dimension
    constexpr double a_p        = 1.792847386e+00;

    /// The charge of proton
    constexpr double z_p        = 1;

    /// The carbon rest mass in GeV
    constexpr double m_c        = 12 * 0.931494027e+00;

    /// The uranium rest mass in GeV
    constexpr double m_u        = 238 * 0.931494027e+00;

    /// The muon rest mass in GeV
    constexpr double m_mu       = 0.10565837;

    /// The h- rest mass in GeV
    constexpr double m_hm       = 0.939277e+00;

    /// The hyodrogen atom rest mass in GeV
    constexpr double m_h        = 1.00782503224 * uam;

    /// The deuteron rest mass in GeV
    constexpr double m_d        = 2*0.931494027e+00;

    /// The xenon rest mass in GeV
    constexpr double m_xe       = 124*0.931494027e+00;

    constexpr double PMASS      = 1.6726231e-27;  // kg

    constexpr double EMASS      = 9.1093897e-31; // kg

    constexpr double PCHARGE    = 1.6021773349e-19; // C

    // alfven current
    constexpr double Ia         = 17.045148e+03;
    // e/mc
    constexpr double e0mc       = 5.86679074042e+02;
    // e/m
    constexpr double e0m        = 1.75881961664e+11;
};

#endif // CLASSIC_Physics_HH