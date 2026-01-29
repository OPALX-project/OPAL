//
// Class Element
//   The base class for all OPAL elements.
//   It implements the common behaviour of elements, it can also be used via
//   dynamic casting to determine whether an object represents an element.
//
//   Each Element object contains a pointer to a CLASSIC beam line element,
//   known as the ``ideal'' element.
//
//   If sharable flag is set, all occurrences of the element are supposed to
//   have the same imperfections.  Thus the assembly is shared when it is used
//   more than once in beam lines or sequences.
//
//   If the sharable flag is not set, each occurrence of the element is supposed
//   to have its own imperfections, but the same ideal representation.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"


Element::~Element()
{}


bool Element::canReplaceBy(Object *object) {
    return (dynamic_cast<Element *>(object) != 0);
}


Element *Element::find(const std::string &name) {
    OpalData *opal = OpalData::getInstance();
    Element *element = dynamic_cast<Element *>(opal->find(name));
    if(element == 0) {
        throw OpalException("Element::find()",
                            "Element \"" + name + "\" not found.");
    }
    return element;
}


const std::string Element::getCategory() const {
    return "ELEMENT";
}


bool Element::shouldTrace() const {
    return false;

}

bool Element::shouldUpdate() const {
    return false;
}


double Element::getEntrance(ReferenceType ref) const {
    switch(ref) {

        case IS_CENTRE:
            return (- getLength() / 2.0);

        case IS_EXIT:
            return (- getLength());

        default:
            return 0.0;
    }
}


double Element::getExit(ReferenceType ref) const {
    switch(ref) {

        case IS_ENTRY:
            return getLength();

        case IS_CENTRE:
            return (getLength() / 2.0);

        default:
            return 0.0;
    }
}


void Element::setShared(bool flag) {
    Object::setShared(flag);
    if(flag) itsClassicElement->makeSharable();
}


Element::Element(int size, const char *name, const char *help):
    Object(size, name, help)
{}


Element::Element(const std::string &name, Element *parent):
    Object(name, parent)
{}
