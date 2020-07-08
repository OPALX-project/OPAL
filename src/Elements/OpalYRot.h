#ifndef OPAL_OpalYRot_HH
#define OPAL_OpalYRot_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalYRot.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// CopSRight: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalYRot
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalYRot
// ------------------------------------------------------------------------
/// The YROT element.

class OpalYRot: public OpalElement {

public:

    /// The attributes of class OpalYRot.
    enum {
        ANGLE = COMMON,  // The rotation angle.
        SIZE             // Total number of attributes.
    };

    /// Exemplar constructor.
    OpalYRot();

    virtual ~OpalYRot();

    /// Make clone.
    virtual OpalYRot *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &);

    /// Update the embedded CLASSIC patch.
    virtual void update();

private:

    // Not implemented.
    OpalYRot(const OpalYRot &);
    void operator=(const OpalYRot &);

    // Clone constructor.
    OpalYRot(const std::string &name, OpalYRot *parent);
};

#endif // OPAL_OpalYRot_HH
