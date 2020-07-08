#ifndef OPAL_OpalRBend3D_HH
#define OPAL_OpalRBend3D_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalRBend3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRBend3D
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:35:09 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalBend.h"

class OpalWake;
class ParticleMatterInteraction;

//
// Class OpalRBend3D
// ------------------------------------------------------------------------
/// The RBEND element.

class OpalRBend3D: public OpalElement {

public:

    enum {
        ANGLE = COMMON,   // The bend angle.
        K0, K0S,          // The multipole coefficients; must be in this order.
        E1,               // The edge angles.
        FMAPFN,           // File name containing on-axis field.
        GAP,              // Full gap of magnet.
        HAPERT,           // Horizontal aperture of magnet.
        DESIGNENERGY,     // the design energy of the particles
        SIZE              // Total number of attributes.
    };

    /// Exemplar constructor.
    OpalRBend3D();

    virtual ~OpalRBend3D();

    /// Make clone.
    virtual OpalRBend3D *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &);

    /// Update the embedded CLASSIC bend.
    virtual void update();

    virtual void print(std::ostream &) const;

private:

    // Not implemented.
    OpalRBend3D(const OpalRBend3D &);
    void operator=(const OpalRBend3D &);

    // Clone constructor.
    OpalRBend3D(const std::string &name, OpalRBend3D *parent);

    OpalWake *owk_m;
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalRBend3D_HH