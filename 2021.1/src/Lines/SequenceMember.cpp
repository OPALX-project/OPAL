// ------------------------------------------------------------------------
// $RCSfile: SequenceMember.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SequenceMember
//   A SequenceMember includes the special data required for OPAL sequences.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:16:16 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Lines/SequenceMember.h"
#include "AbsBeamline/ElementBase.h"


// Class SequenceMember
// ------------------------------------------------------------------------

SequenceMember::SequenceMember():
    FlaggedElmPtr(),
    itsPosition(0.0), itsFlag(ABSOLUTE), itsType(UNKNOWN), OpalElement()
{}


SequenceMember::SequenceMember(const SequenceMember &rhs):
    FlaggedElmPtr(rhs),
    itsPosition(rhs.itsPosition),
    itsFlag(rhs.itsFlag),
    itsType(rhs.itsType),
    OpalElement(rhs.OpalElement)
{}


SequenceMember::~SequenceMember()
{}


void SequenceMember::setLength(double drift) {
    if(itsType == GENERATED) getElement()->setElementLength(drift);
}
