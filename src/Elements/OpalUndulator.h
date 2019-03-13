#ifndef OPAL_OpalUndulator_HH
#define OPAL_OpalUndulator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalUndulator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalUndulator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class BoundaryGeometry;

// Class OpalUndulator
// ------------------------------------------------------------------------
/// The DRIFT element.

class OpalWake;
class ParticleMatterInteraction;

class OpalUndulator: public OpalElement {

public:

    enum {
         GEOMETRY = COMMON,       // geometry of boundary, one more enum member besides the common ones in OpalElement.
		 NSLICES,	  // The number of slices / steps per element for map tracking
	 K,                       // The undulator parameter
	 FNAME,
	 SIZE

    };
    /// Exemplar constructor.
    OpalUndulator();

    virtual ~OpalUndulator();

    /// Make clone.
    virtual OpalUndulator *clone(const std::string &name);

    /// Test for drift.
    //  Return true.
    virtual bool isUndulator() const;

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:

    // Not implemented.
    OpalUndulator(const OpalUndulator &);
    void operator=(const OpalUndulator &);

    // Clone constructor.
    OpalUndulator(const std::string &name, OpalUndulator *parent);

    OpalWake *owk_m;
    ParticleMatterInteraction *parmatint_m;
    BoundaryGeometry *obgeo_m;
};

#endif // OPAL_OpalUndulator_HH
