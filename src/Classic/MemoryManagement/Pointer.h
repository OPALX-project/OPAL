#ifndef CLASSIC_Pointer_HH
#define CLASSIC_Pointer_HH

// ------------------------------------------------------------------------
// $RCSfile: Pointer.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Pointer<Object>
//
// ------------------------------------------------------------------------
// Class category: MemoryManagement
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Template class Pointer
// ------------------------------------------------------------------------
/// Reference-counted pointer.
//  This template class implements a pointer with reference counting in
//  collaboration with the class RCObject. It is modelled after the class
//  RCPtr described in:
//  {center}
//  Scott Meyers, More Effective C++, Addison Wesley, 1996, pg. 195.
//  {/center}
//  When the pointer goes out of scope, the pointee's reference count is
//  decremented, and, if it becomes zero, the pointee is deleted.
//  When a Pointer is assigned, the pointee's reference count is incremented.

template < class Object >
class Pointer {

public:

    /// Default constructor.
    //  Set the pointer to nullptr.
    Pointer();

    /// Copy constructor.
    //  Increment reference count of pointee.
    Pointer(const Pointer &);

    /// Constructor.
    //  Grab an object just created by "new" and set its reference count to one.
    Pointer(Object *);

    /// Destructor.
    //  Decrement reference count and delete object if it reaches zero.
    ~Pointer();

    /// Assign.
    //  Increment reference count of pointee.
    Pointer &operator=(const Pointer &);

    /// Assign.
    //  Grab an object just created by "new" and set its reference count to one.
    Pointer &operator=(Object *);

    /// Delegation operator.
    //  Return a C-type pointer to the object.
    Object *operator->() const;

    /// Dereferencing operator.
    //  Return a reference to the object.
    Object &operator*() const;

    /// Pointer equality.
    bool operator==(const Pointer &) const;

    /// Pointer inequality.
    bool operator!=(const Pointer &) const;

    /// Test for validity.
    //  Return true, if the pointer is not nullptr.
    bool isValid() const;

    /// Force unique.
    //  Copy the object and make *this point to it.
    void unique();

private:

    // The object pointed at.
    Object *object;
};


// Inline functions for template class Pointer
// ------------------------------------------------------------------------

template < class Object >
inline Pointer<Object>::Pointer():
    object(0)
{}


template < class Object >
inline Pointer<Object>::Pointer(const Pointer<Object> &rhs):
    object(rhs.object) {
    if(object) object->addReference();
}


template < class Object >
inline Pointer<Object>::Pointer(Object *obj):
    object(obj) {
    if(object) object->addReference();
}


template < class Object >
inline Pointer<Object>::~Pointer() {
    if(object  &&  object->removeReference() <= 0) delete object;
}


template < class Object >
inline Pointer<Object> &Pointer<Object>::operator=(const Pointer &rhs) {
    if(object != rhs.object) {
        if(object != 0  &&  object->removeReference() <= 0) delete object;
        object = rhs.object;
        if(object) object->addReference();
    }

    return *this;
}


template < class Object >
inline Pointer<Object> &Pointer<Object>::operator=(Object *obj) {
    if(object != obj) {
        if(object  &&  object->removeReference() <= 0) delete object;
        object = obj;
        if(object) object->addReference();
    }
    return *this;
}


template < class Object >
inline Object *Pointer<Object>::operator->() const {
    return object;
}


template < class Object >
inline Object &Pointer<Object>::operator*() const {
    return *object;
}


template < class Object >
inline bool Pointer<Object>::operator==(const Pointer &rhs) const {
    return object == rhs.object;
}


template < class Object >
inline bool Pointer<Object>::operator!=(const Pointer &rhs) const {
    return object != rhs.object;
}


template < class Object >
inline bool Pointer<Object>::isValid() const {
    return object != 0;
}


template < class Object >
inline void Pointer<Object>::unique() {
    if(object != 0  &&  object->isShared()) {
        // The static cast is safe, since the declaration says
        // "Object *object;"
        object->removeReference();
        object = static_cast<Object *>(object->clone());
        object->addReference();
    }
}

#endif // CLASSIC_Pointer_HH
