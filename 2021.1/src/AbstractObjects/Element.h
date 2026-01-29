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
#ifndef OPAL_Element_HH
#define OPAL_Element_HH

#include "AbstractObjects/Object.h"
#include "AbsBeamline/ElementBase.h"
#include "MemoryManagement/Pointer.h"


class Element: public Object {

public:

    /// Reference for element positioning.
    //  Used in the SEQUENCE command.
    enum ReferenceType {
        IS_ENTRY,         // Reference point is at element entrance.
        IS_CENTRE,        // Reference point is at element centre
        IS_EXIT           // Reference point is at element exit.
    };

    virtual ~Element();

    /// Test if replacement is allowed.
    //  Return true, if the replacement is also an Element.
    virtual bool canReplaceBy(Object *object);

    /// Find named Element.
    //  If an element with the name [b]name[/b] exists,
    //  return a pointer to that element.
    //  If no such element exists, throw [b]OpalException[/b].
    static Element *find(const std::string &name);

    /// Return the object category as a string.
    //  Return the string "ELEMENT".
    virtual const std::string getCategory() const;

    /// Trace flag.
    //  If true, the object's execute() function should be traced.
    //  Always false for elements.
    virtual bool shouldTrace() const;

    /// Update flag.
    //  If true, the data structure should be updated before calling execute().
    //  Always false for elements.
    virtual bool shouldUpdate() const;


    /// Return element length.
    virtual double getLength() const = 0;

    /// Return arc length from origin to entrance (negative !).
    virtual double getEntrance(ReferenceType) const;

    /// Return arc length from origin to exit (positive !).
    virtual double getExit(ReferenceType) const;

    /// Set shared flag.
    //  If true, all references to this name are to the same object.
    virtual void setShared(bool);

    /// Return the embedded CLASSIC element.
    //  Return a pointer to the embedded CLASSIC ElementBase
    inline ElementBase *getElement() const;

    /// Assign new CLASSIC element.
    inline void setElement(ElementBase *);

protected:

    /// Constructor for exemplars.
    Element(int size, const char *name, const char *help);

    /// Constructor for clones.
    Element(const std::string &name, Element *parent);

private:

    // Not implemented.
    Element();
    Element(const Element &);
    void operator=(const Element &);

    // The embedded CLASSIC element.
    Pointer<ElementBase> itsClassicElement;
};


// Inline functions.
// ------------------------------------------------------------------------

inline ElementBase *Element::getElement() const {
    return &*itsClassicElement;
}


inline void Element::setElement(ElementBase *base) {
    itsClassicElement = base;
}

#endif // OPAL_Element_HH
