#ifndef OPAL_OpalSRot_HH
#define OPAL_OpalSRot_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSRot.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// CopSRight: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSRot
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalSRot
// ------------------------------------------------------------------------
/// The SROT element.

class OpalSRot: public OpalElement {

public:

    /// The attributes of class OpalSRot.
    enum {
        ANGLE = COMMON,  // The rotation angle.
        SIZE             // Total number of attributes.
    };

    /// Exemplar constructor.
    OpalSRot();

    virtual ~OpalSRot();

    /// Make clone.
    virtual OpalSRot *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &);

    /// Update the embedded CLASSIC patch.
    virtual void update();

private:

    // Not implemented.
    OpalSRot(const OpalSRot &);
    void operator=(const OpalSRot &);

    // Clone constructor.
    OpalSRot(const std::string &name, OpalSRot *parent);
};

#endif // OPAL_OpalSRot_HH
