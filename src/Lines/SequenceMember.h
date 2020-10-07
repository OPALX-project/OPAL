#ifndef OPAL_SequenceMember_HH
#define OPAL_SequenceMember_HH

// ------------------------------------------------------------------------
// $RCSfile: SequenceMember.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SequenceMember
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:16:16 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Beamlines/FlaggedElmPtr.h"
#include "MemoryManagement/Pointer.h"
#include "AbstractObjects/Element.h"

// Class SequenceMember
// ------------------------------------------------------------------------
/// A member of a SEQUENCE.

class SequenceMember: public FlaggedElmPtr {

public:

    // The code for different SequenceMember types.
    enum MemberType {
        UNKNOWN   = 0,
        GLOBAL    = 1,
        LOCAL     = 2,
        GENERATED = 3
    };

    // The code for different position types.
    enum PositionType {
        ABSOLUTE  = 0,
        FROM      = 1,
        PREVIOUS  = 2,
        NEXT      = 3,
        BEGIN     = 4,
        END       = 5,
        DRIFT     = 6,
        IMMEDIATE = 7
    };

    SequenceMember();
    SequenceMember(const SequenceMember &);
    ~SequenceMember();

    /// Store the drift length for a generated drift.
    void setLength(double drift);


    /// The position attribute ("AT" or "DRIFT").
    double itsPosition;

    /// Flag word.
    //  Specific to the sequence parser, tells how the element is related
    //  to the preceding one.
    PositionType itsFlag;

    /// Type word.
    //  Tells how the element is defined.
    MemberType itsType;

    // ada 4.5 2000 to speed up matching, add a pointer to
    // opal elements in order to avoid searching the opal elements
    Pointer<Element> OpalElement;


private:

    // Not implemented.
    void operator=(const SequenceMember &);
};

#endif // OPAL_SequenceMember_HH
