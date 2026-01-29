#ifndef OPAL_OpalMarker_HH
#define OPAL_OpalMarker_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalMarker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMarker
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalMarker
// ------------------------------------------------------------------------
/// The MARKER element.

class OpalMarker: public OpalElement {

public:

    /// Exemplar constructor.
    OpalMarker();

    virtual ~OpalMarker();

    /// Make clone.
    virtual OpalMarker *clone(const std::string &name);

    /// Print the element.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Update the embedded CLASSIC marker.
    virtual void update();

private:

    // Not implemented.
    OpalMarker(const OpalMarker &);
    void operator=(const OpalMarker &);

    // Clone constructor.
    OpalMarker(const std::string &name, OpalMarker *parent);
};

#endif // OPAL_OpalMarker_HH
