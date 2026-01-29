#ifndef CLASSIC_FlaggedElmPtr_HH
#define CLASSIC_FlaggedElmPtr_HH

// ------------------------------------------------------------------------
// $RCSfile: FlaggedElmPtr.h,v $
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

#include "Beamlines/ElmPtr.h"


// Class FlaggedElmPtr
// ------------------------------------------------------------------------
/// A section of a beam line.
//  A FlaggedBeamline is built as a list of FlaggedElmPtr objects.
//  A FlaggedElmPtr contains two flags:
//  [dl]
//  [dt]reflected [dd]If true, this section of the line is reversed.
//  [dt]selected  [dd]Can be set to indicate a selected status for algorithms.
//  [/dl]

class FlaggedElmPtr: public ElmPtr {

public:

    /// Constructor.
    //  Assign an element to this beamline position.
    explicit FlaggedElmPtr
    (const ElmPtr &, bool reflected = false, bool selected = false);

    /// Copy constructor.
    FlaggedElmPtr(const FlaggedElmPtr &);

    FlaggedElmPtr();
    virtual ~FlaggedElmPtr();

    /// Apply visitor.
    //  If any error occurs, this method throws an exception.
    virtual void accept(BeamlineVisitor &) const;

    /// Get clone counter.
    //  Return the value set by the last call to [tt]setCounter()[/tt].
    //  See [b]FlaggedElmPtr::itsCounter[/b] for details.
    inline int getCounter() const;

    /// Get reflection flag.
    inline bool getReflectionFlag() const;

    /// Get selection flag.
    inline bool getSelectionFlag() const;

    /// Set clone counter.
    //  The value set is the entire responsibility of any algorithm using it.
    //  See [b]FlaggedElmPtr::itsCounter[/b] for details.
    inline void setCounter(int) const;

    /// Set reflection flag.
    inline void setReflectionFlag(bool flag) const;

    /// Get selection flag.
    inline void setSelectionFlag(bool flag) const;

protected:

    /// Clone counter.
    //  This value can be set and interrogated by an algorithm.
    //  It is not used by the CLASSIC library.
    mutable int itsCounter;

    /// The reflection flag.
    //  If true, the portion of the line pointed at by [b]itsElement[/b]
    //  is reflected, i.e. its elements occur in reverse order.
    mutable bool isReflected;

    /// The selection flag.
    //  This flag can be set to indicate a ``selected'' status for certain
    //  algorithms.
    mutable bool isSelected;
};


inline int FlaggedElmPtr::getCounter() const {
    return itsCounter;
}


inline bool FlaggedElmPtr::getReflectionFlag() const {
    return isReflected;
}


inline bool FlaggedElmPtr::getSelectionFlag() const {
    return isSelected;
}


inline void FlaggedElmPtr::setCounter(int count) const {
    itsCounter = count;
}


inline void FlaggedElmPtr::setReflectionFlag(bool flag) const {
    isReflected = flag;
}


inline void FlaggedElmPtr::setSelectionFlag(bool flag) const {
    isSelected = flag;
}

#endif // CLASSIC_FlaggedElmPtr_HH
