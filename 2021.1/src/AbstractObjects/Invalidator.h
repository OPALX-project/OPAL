#ifndef OPAL_Invalidator_HH
#define OPAL_Invalidator_HH

// ------------------------------------------------------------------------
// $RCSfile: Invalidator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Invalidator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:40:35 $
// $Author: opal $
//
// ------------------------------------------------------------------------


// Class Invalidator
// ------------------------------------------------------------------------
/// Abstract base class for references which must be invalidated when an
//  object goes out of scope.

class Invalidator {

public:

    /// Force re-evaluation.
    //  Set an internal flag so as to force re-evaluation of any expression
    //  when the value is referred next time.
    virtual void invalidate();
};

#endif // OPAL_Invalidator_HH
