#ifndef OPAL_OpalQuadrupole_HH
#define OPAL_OpalQuadrupole_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalQuadrupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalQuadrupole
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

// Class OpalQuadrupole
// ------------------------------------------------------------------------
/// The QUADRUPOLE element.

class OpalQuadrupole: public OpalElement {

public:

    /// The attributes of class OpalQuadrupole.
    enum {
        K1 = COMMON,  // The normal quadrupole coefficient.
        DK1,          // The normal quadupole coefficient error.
        K1S,          // The skew quadrupole coefficient.
        DK1S,         // The skew quadrupole coefficient error.
        NSLICES,      // The number of slices / steps per element for map tracking
        SIZE
    };

    /// Exemplar constructor.
    OpalQuadrupole();

    virtual ~OpalQuadrupole();

    /// Make clone.
    virtual OpalQuadrupole *clone(const std::string &name);

    /// Print the quadrupole.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalQuadrupole(const OpalQuadrupole &);
    void operator=(const OpalQuadrupole &);

    // Clone constructor.
    OpalQuadrupole(const std::string &name, OpalQuadrupole *parent);

    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalQuadrupole_HH