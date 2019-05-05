#ifndef OPAL_OpalCavity_HH
#define OPAL_OpalCavity_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalCavity.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCavity
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

// Class OpalCavity
// ------------------------------------------------------------------------
/// The RFCAVITY element.

class OpalWake;
class BoundaryGeometry;

class OpalCavity: public OpalElement {

public:

    /// The attributes of class OpalCavity.
    enum {
        VOLT = COMMON,  // The peak voltage.
        DVOLT,          // The peak voltage error.
        GEOMETRY,       // geometry of boundary
        FREQ,           // The RF frequency.
        LAG,            // The phase lag.
        DLAG,           // The phase lag error.
        HARMON,         // The harmonic number.
        BETARF,         // The beta_RF.
        PG,             // The RF power.
        ZSHUNT,         // The shunt impedance.
        TFILL,          // The filling time.
        FMAPFN,         // The filename of the fieldmap
        FAST,           // Faster but less accurate
        APVETO,         // Do not use this cavity in the Autophase procedure
        RMIN,           // Minimal Radius
        RMAX,           // Maximal Radius
        ANGLE,          // the azimuth position of the cavity
        PDIS,           // perpendicular distance from symmetric line of cavity gap to machine center
        GAPWIDTH,       // constant gap width of cavity
        PHI0,           // initial phase of cavity
        DESIGNENERGY,   // The mean kinetic energy at exit
        PHASE_MODEL,    // time dependent phase
        AMPLITUDE_MODEL,// time dependent amplitude
        FREQUENCY_MODEL,// time dependent frequency
        SIZE
    };

    /// Exemplar constructor.
    OpalCavity();

    virtual ~OpalCavity();

    /// Make clone.
    virtual OpalCavity *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalCavity(const OpalCavity &);
    void operator=(const OpalCavity &);

    // Clone constructor.
    OpalCavity(const std::string &name, OpalCavity *parent);

    OpalWake *owk_m;

    BoundaryGeometry *obgeo_m;


};

#endif // OPAL_OpalCavity_HH