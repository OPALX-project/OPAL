// ------------------------------------------------------------------------
// $RCSfile: Replacer.cpp,v $
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

#include "Lines/Replacer.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"

class Element;


// Class Replacer
// ------------------------------------------------------------------------

Replacer::Replacer(const Beamline &beamline,
                   const std::string &name,
                   ElementBase *elm):
    DefaultVisitor(beamline, false, false),
    itsName(name),
    newBase(elm)
{}


Replacer::~Replacer()
{}


void Replacer::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    // Find proper OPAL element.
    const std::string &name = fep.getElement()->getName();

    // Do the required operations on the beamline or element.
    if(name == itsName) {
        ElementBase *base = newBase->copyStructure();
        const_cast<FlaggedElmPtr &>(fep).setElement(base);
    } else {
        DefaultVisitor::visitFlaggedElmPtr(fep);
    }
}
