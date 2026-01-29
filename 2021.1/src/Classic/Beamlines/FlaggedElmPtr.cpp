// ------------------------------------------------------------------------
// $RCSfile: FlaggedElmPtr.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FlaggedElmPtr
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/FlaggedElmPtr.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/LogicalError.h"
#include <stdexcept>


// Class FlaggedElmPtr
// ------------------------------------------------------------------------

FlaggedElmPtr::FlaggedElmPtr
(const ElmPtr &elem, bool reflected, bool selected):
    ElmPtr(elem),
    itsCounter(0),
    isReflected(reflected),
    isSelected(selected)
{}


FlaggedElmPtr::FlaggedElmPtr():
    ElmPtr(),
    itsCounter(0),
    isReflected(false),
    isSelected(false)
{}


FlaggedElmPtr::FlaggedElmPtr(const FlaggedElmPtr &rhs):
    ElmPtr(rhs),
    itsCounter(rhs.itsCounter),
    isReflected(rhs.isReflected),
    isSelected(rhs.isSelected)
{}


FlaggedElmPtr::~FlaggedElmPtr()
{}


void FlaggedElmPtr::accept(BeamlineVisitor &v) const {
    v.visitFlaggedElmPtr(*this);
}
