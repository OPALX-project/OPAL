#ifndef CLASSIC_RCObject_HH
#define CLASSIC_RCObject_HH

// ------------------------------------------------------------------------
// $RCSfile: RCObject.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RCObject
//
// ------------------------------------------------------------------------
// Class category: MemoryManagement
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class RCObject
// ------------------------------------------------------------------------
/// Abstract base class for reference counted objects.
//  It collaborates with the templace class Pointer<Object>.
//  It is modelled after the class RCObject described in:
//  {center}
//  Scott Meyers, More Effective C++, Addison Wesley, 1996, pg. 195.
//  {/center}
//  The idiom [tt]delete this[/tt] has been avoided by moving the deletion
//  to class Pointer<Object>.
//  All derived classes must implement a clone() method, to allow the
//  Pointer class to copy the object pointed at.
//  All constructors, the destructor, and the assignment operators are
//  protected, since a stand-alone RCObject makes no sense.


class RCObject {

public:

    /// Increment reference count.
    //  Return the new value of the reference count.
    int addReference() const;

    /// Decrement the reference count.
    //  Return the new value of the reference count.
    int removeReference() const;

    /// Test for sharing.
    //  Return true, if the pointee has more than one reference.
    bool isShared() const;


protected:

    /// Default constructor.
    //  By default a new object is sharable.
    RCObject();

    /// Copy constructor.
    //  The copy inherits the sharable state.
    RCObject(const RCObject &);

    // Destructor is pure, but implemented, to make class abstract.
    virtual ~RCObject() = 0;

    RCObject &operator=(const RCObject &right);

private:

    // The object's reference count.
    // The value is mutable, since it is not really part of the objects state.
    mutable int refCount;
};


// Inline functions for class RCObject
// ------------------------------------------------------------------------

inline int RCObject::addReference() const {
    return ++refCount;
}


inline int RCObject::removeReference() const {
    return --refCount;
}


inline bool RCObject::isShared() const {
    return refCount > 1;
}


inline RCObject::RCObject():
    refCount(0)
{}


inline RCObject::RCObject(const RCObject &/*rhs*/):
    refCount(0)
{}


inline RCObject::~RCObject()
{}


inline RCObject &RCObject::operator=(const RCObject &/*rhs*/) {
    return *this;
}

#endif // CLASSIC_RCObject_HH
