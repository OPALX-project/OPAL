#ifndef MAD_Flagger_HH
#define MAD_Flagger_HH

// ------------------------------------------------------------------------
// $RCSfile: Flagger.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Flagger
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"

class Beamline;


// Class Flagger
// ------------------------------------------------------------------------
/// Set/reset all selection flags in a beam line built from FlaggedElmPtr
//  objects.

class Flagger: public DefaultVisitor {

public:

    /// Constructor.
    //  Attach this visitor to [b]bl[/b], remember the [b]set[/b] flag.
    Flagger(const Beamline &bl, bool set);

    virtual ~Flagger();

    /// Set selection flag in the given FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

private:

    // Not implemented.
    Flagger();
    Flagger(const Flagger &);
    void operator=(const Flagger &);

    // The flag value to be applied.
    bool flag;
};

#endif // MAD_Flagger_HH
