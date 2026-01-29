#ifndef OPAL_Replacer_HH
#define OPAL_Replacer_HH

// ------------------------------------------------------------------------
// $RCSfile: Replacer.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Replacer
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include <string>

class Beamline;
class ElementBase;


// Class Replacer
// ------------------------------------------------------------------------
/// Replace all references to named element by a new version.

class Replacer: public DefaultVisitor {

public:

    /// Constructor.
    //  Attach visitor to a beamline, remember the range.
    Replacer(const Beamline &, const std::string &name, ElementBase *elm);

    virtual ~Replacer();

    /// Apply the visitor to an FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

private:

    // Not implemented.
    Replacer();
    Replacer(const Replacer &);
    void operator=(const Replacer &);

    // The name to be replaced.
    const std::string itsName;

    // The new element replacing the old one.
    ElementBase *newBase;
};

#endif // OPAL_Replacer_HH
