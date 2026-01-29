#ifndef OPAL_RangeRep_HH
#define OPAL_RangeRep_HH

// ------------------------------------------------------------------------
// $RCSfile: RangeRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RangeRep
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/PlaceRep.h"

class FlaggedElmPtr;


// Class RangeRep
// ------------------------------------------------------------------------
/// Representation of a range within a beam line or sequence.
//  This representation holds two places, defining the beginning and ending
//  of the range.  It keeps track wether we are before, within, or after
//  the range.
//  The value may also be "FULL", denoting the complete beam line.

class RangeRep {

public:

    /// Default constructor.
    //  Construct range for full line (value "FULL").
    RangeRep();

    /// Constructor from two given places.
    RangeRep(PlaceRep &f, PlaceRep &l);

    RangeRep(const RangeRep &);
    ~RangeRep();
    const RangeRep &operator=(const RangeRep &);

    /// Initialise data for search.
    //  Sets the internal state to the beginning of the beam line.
    void initialize();

    /// Test for active range.
    //  Return true, if we are within the specified range.
    bool isActive() const;

    /// Enter an element or line.
    //  Sets the internal state to active, when we enter the specified range.
    void enter(const FlaggedElmPtr &) const;

    /// Leave an element or line.
    //  Sets the internal state to inactive, when we leave the specified range.
    void leave(const FlaggedElmPtr &) const;

    /// Print in input format.
    void print(std::ostream &os) const;

private:

    // The begin/end Place structures.
    PlaceRep first;
    PlaceRep last;

    // Flag for full range. If true the places are not looked at.
    bool fullRange;

    // If true we are within this range.
    mutable bool status;
};


inline std::ostream &operator<<(std::ostream &os, const RangeRep &val) {
    val.print(os);
    return os;
}

#endif // OPAL_RangeRep_HH
