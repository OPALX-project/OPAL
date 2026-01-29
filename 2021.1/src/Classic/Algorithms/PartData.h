#ifndef MAD_PartData_HH
#define MAD_PartData_HH

// ------------------------------------------------------------------------
// $RCSfile: PartData.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

// Class PartData
// ------------------------------------------------------------------------
/// Particle reference data.
//  This class encapsulates the reference data for a beam:
//  [UL]
//  [LI]charge per particle expressed in proton charges,
//  [LI]mass per particle expressed in eV,
//  [LI]reference momentum per particle expressed in eV.
//  [/UL]
//  The copy constructor, destructor, and assignment operator generated
//  by the compiler perform the correct operation.  For speed reasons
//  they are not implemented.

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


    /// Set reference mass expressed in eV/c^2
    inline void setM(double m){mass = m;}

    /// Set reference charge expressed in proton charges,  
    inline void setQ(double q) {charge = q;}

protected:

    // The reference information.
    double charge;   // Particle charge.
    double mass;     // Particle mass.
    double beta;     // particle velocity divided by c.
    double gamma;    // particle energy divided by particle mass
};


// Inline functions.
// ------------------------------------------------------------------------

inline double PartData::getQ() const {
    return charge;
}


inline double PartData::getM() const {
    return mass;
}


inline double PartData::getP() const {
    return beta * gamma * mass;
}


inline double PartData::getE() const {
    return gamma * mass;
}


inline double PartData::getBeta() const {
    return beta;
}


inline double PartData::getGamma() const {
    return gamma;
}

#endif // MAD_PartData_HH
