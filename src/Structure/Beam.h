#ifndef OPAL_Beam_HH
#define OPAL_Beam_HH

// ------------------------------------------------------------------------
// $RCSfile: Beam.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Beam
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include <ostream>

class Inform;

// Class Beam
// ------------------------------------------------------------------------
/// The BEAM definition.
//  A BEAM definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class Beam: public Definition {

public:

    /// Exemplar constructor.
    Beam();

    virtual ~Beam();

    /// Test if replacement is allowed.
    //  Can replace only by another BEAM.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual Beam *clone(const std::string &name);

    /// Check the BEAM data.
    virtual void execute();

    /// Find named BEAM.
    static Beam *find(const std::string &name);

    /// Return emittance for mode 1.
    double getEX() const;

    /// Return emittance for mode 2.
    double getEY() const;

    /// Return emittance for mode 3.
    double getET() const;

    //ff => get gamma value
    double getGamma() const;

    //ff => get PC value
    double getPC() const;

    /// Return the number of (macro)particles
    size_t getNumberOfParticles();

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;

    /// Return the beam current in A
    double getCurrent() const;

    /// Return the charge number in elementary charge
    double getCharge() const;

    /// Return the beam frequency in MHz
    double getFrequency() const;

    /// Return Particle's name
    std::string getParticleName() const;

    /// Return Particle's rest mass in GeV
    double getMass() const;

    /// Store emittance for mode 1.
    void setEX(double);

    /// Store emittance for mode 2.
    void setEY(double);

    /// Store emittance for mode 3.
    void setET(double);

    /// Update the BEAM data.
    virtual void update();

    void print(std::ostream &os) const;

private:

    // Not implemented.
    Beam(const Beam &);
    void operator=(const Beam &);

    // Clone constructor.
    Beam(const std::string &name, Beam *parent);

    // The particle reference data.
    PartData reference;

    // The converstion from GeV to eV.
    static const double energy_scale;
};

inline std::ostream &operator<<(std::ostream &os, const Beam &b) {
    b.print(os);
    return os;
}


#endif // OPAL_Beam_HH
