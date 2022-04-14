#ifndef OPAL_OpalOutputPlane_H
#define OPAL_OpalOutputPlane_H

#include "Elements/OpalElement.h"

// Class OpalOutputPlane
// ------------------------------------------------------------------------
/// Interface for output plane element.


class OpalOutputPlane: public OpalElement {

public:

    /// The attributes of class OpalOutputPlane.
    enum {
        CENTRE_X = COMMON, // x coordinate of plane centre
        CENTRE_Y,          // y coordinate of plane centre
        CENTRE_Z,          // z coordinate of plane centre
        NORMAL_X,          // x coordinate of plane normal
        NORMAL_Y,          // y coordinate of plane normal
        NORMAL_Z,          // z coordinate of plane normal
        TOLERANCE,         // tolerance on position estimate
        REFERENCE_ALIGNMENT_PARTICLE, // if true, centres on reference particle
        WIDTH,             // maximum allowed horizontal displacement from centre
        HEIGHT,            // maximum allowed vertical displacement from centre
        RADIUS,            // maximum allowed radial displacement from centre
        ALGORITHM,         // algorithm used to find crossing-point
        MATERIAL,          // assume a material for energy loss and scattering during crossing
        MATERIAL_THICKNESS, // algorithm used to calculate material thickness
        SIZE
    };
    /// Exemplar constructor.
    OpalOutputPlane();

    virtual ~OpalOutputPlane();

    /// Make clone.
    virtual OpalOutputPlane *clone(const std::string &name);

    /// Fill in all registered attributes.
    //virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalOutputPlane(const OpalOutputPlane &);
    void operator=(const OpalOutputPlane &);

    // Clone constructor.
    OpalOutputPlane(const std::string &name, OpalOutputPlane *parent);
};

#endif // OPAL_OpalOutputPlane_H
