#ifndef OPAL_OpalBeamStripping_HH
#define OPAL_OpalBeamStripping_HH

// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Class: OpalBeamStripping
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;
// Class OpalBeamStripping
// ------------------------------------------------------------------------
/// The BEAMSTRIPPING element.

class OpalBeamStripping: public OpalElement {

public:

    /// The attributes of class OpalBeamStripping.
    enum {
        PRESSURE = COMMON,   // Pressure in mbar
        TEMPERATURE,         // Temperature of residual gas
        PMAPFN,              // The filename of the mid-plane pressure map
        PSCALE,              // A scalar to scale the P-field
        GAS,                 // Type of gas for residual vaccum
        STOP,                // whether the secondary particles are tracked
        SIZE
    };

    /// Exemplar constructor.
    OpalBeamStripping();

    virtual ~OpalBeamStripping();

    /// Make clone.
    virtual OpalBeamStripping *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC beam stripping.
    virtual void update();

private:

    // Not implemented.
    OpalBeamStripping(const OpalBeamStripping &);
    void operator=(const OpalBeamStripping &);

    // Clone constructor.
    OpalBeamStripping(const std::string &name, OpalBeamStripping *parent);
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalBeamStripping_HH