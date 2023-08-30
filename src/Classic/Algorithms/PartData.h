//
// Class PartData
//   PartData represents a set of reference values for use in algorithms.
//
// Copyright (c) 200x - 2023, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef MAD_PartData_HH
#define MAD_PartData_HH

/** Class PartData
 * ------------------------------------------------------------------------
 *  Particle reference data.
 *  This class encapsulates the reference data for a beam:
 *  [UL]
 *  [LI]charge per particle expressed in proton charges,
 *  [LI]mass per particle expressed in eV,
 *  [LI]reference momentum per particle expressed in eV.
 *  [LI]momentumTolerance Fractional tolerance to deviations in the distribution
 *                        compared to the reference data at initialisation
 *                        If negative, no tolerance checking is done.
 *  [/UL]
 *  The copy constructor, destructor, and assignment operator generated
 *  by the compiler perform the correct operation.  For speed reasons
 *  they are not implemented.
 */

class PartData {

public:
    /// Constructor.
    //  Inputs are:
    //  [DL]
    //  [DT]charge[DD]The charge per particle in proton charges.
    //  [DT]mass[DD]The particle mass in eV.
    //  [DT]momentum[DD]The reference momentum per particle in eV.
    //  [/DL]
    PartData(double charge, double mass, double momentum);

    PartData();

    /// The constant charge per particle.
    double getQ() const;

    /// The constant mass per particle.
    double getM() const;

    /// The constant reference momentum per particle.
    double getP() const;

    /// The constant reference Energy per particle.
    double getE() const;

    /// The relativistic beta per particle.
    double getBeta() const;

    /// The relativistic gamma per particle.
    double getGamma() const;

    /// Get the momentum tolerance
    double getMomentumTolerance() const;

    /// Set reference charge expressed in proton charges.
    void setQ(double q);

    /// Set reference mass expressed in eV/c^2.
    void setM(double m);

    /// Set reference momentum.
    //  Input is the momentum in eV.
    void setP(double p);

    /// Set reference energy.
    //  Input is the energy in eV.
    void setE(double E);

    /// Set beta.
    //  Input is the relativistic beta = v/c.
    void setBeta(double beta);

    /// Set gamma.
    //  Input is the relativistic gamma = E/(m*c*c).
    void setGamma(double gamma);

    /// Set the momentum tolerance.
    void setMomentumTolerance(double tolerance);

protected:
    // The reference information.
    double charge_m;   // Particle charge.
    double mass_m;     // Particle mass.
    double beta_m;     // Particle velocity divided by c.
    double gamma_m;    // Particle energy divided by particle mass.
    double momentumTolerance_m = 1e-2; // Tolerance to momentum deviations.
};


// Inline functions.
// ------------------------------------------------------------------------
inline void PartData::setM(double m) {
    mass_m = m;
}

inline void PartData::setQ(double q) {
    charge_m = q;
}

inline double PartData::getQ() const {
    return charge_m;
}

inline double PartData::getM() const {
    return mass_m;
}

inline double PartData::getP() const {
    return beta_m * gamma_m * mass_m;
}

inline double PartData::getE() const {
    return gamma_m * mass_m;
}

inline double PartData::getBeta() const {
    return beta_m;
}

inline double PartData::getGamma() const {
    return gamma_m;
}

inline double PartData::getMomentumTolerance() const {
    return momentumTolerance_m;
}

#endif // MAD_PartData_HH
