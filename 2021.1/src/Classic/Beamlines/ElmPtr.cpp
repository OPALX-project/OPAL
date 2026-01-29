// ------------------------------------------------------------------------
// $RCSfile: ElmPtr.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElmPtr
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/ElmPtr.h"


// Typedef ElmPtr
// ------------------------------------------------------------------------

ElmPtr::ElmPtr():
    itsElement()
{}


ElmPtr::ElmPtr(const ElmPtr &rhs):
    itsElement(rhs.itsElement)
{}


ElmPtr::ElmPtr(ElementBase *elem):
    itsElement(elem)
{}


ElmPtr::~ElmPtr()
{}


void ElmPtr::accept(BeamlineVisitor &visitor) const {
    itsElement->accept(visitor);
}
