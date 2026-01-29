#ifndef CLASSIC_ElmPtr_HH
#define CLASSIC_ElmPtr_HH

// ------------------------------------------------------------------------
// $RCSfile: ElmPtr.h,v $
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
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"
#include "MemoryManagement/Pointer.h"


// Class ElmPtr.
// ------------------------------------------------------------------------
/// A section of a beam line.
//  A beam line is built as a list of ElmPtr.

class ElmPtr {

public:

    ElmPtr();
    ElmPtr(const ElmPtr &);
    ElmPtr(ElementBase *);
    virtual ~ElmPtr();

    /// Apply visitor.
    //  If any error occurs, this method throws an exception.
    virtual void accept(BeamlineVisitor &) const;

    /// Get the element pointer.
    inline ElementBase *getElement() const;

    /// Set the element pointer.
    inline void setElement(ElementBase *);

protected:

    // The pointer to the element.
    Pointer<ElementBase> itsElement;
};


inline ElementBase *ElmPtr::getElement() const {
    return &*itsElement;
}


inline void ElmPtr::setElement(ElementBase *elem) {
    itsElement = elem;
}

#endif // CLASSIC_ElmPtr_HH
