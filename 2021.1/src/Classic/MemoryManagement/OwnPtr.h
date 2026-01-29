#ifndef CLASSIC_OwnPtr_HH
#define CLASSIC_OwnPtr_HH

// ------------------------------------------------------------------------
// $RCSfile: OwnPtr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: OwnPtr<Object>
//
// ------------------------------------------------------------------------
// Class category: MemoryManagement
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Template class OwnPtr
// ------------------------------------------------------------------------
/// A pointer which owns the object pointed at.
//  When the pointer goes out of scope, the pointee is deleted.
//  When an OwnPtr is assigned the ownership is transferred to the target.

template <class Object>
class OwnPtr {

public:

    /// Default constructor.
    //  Set the pointer to NULL.
    OwnPtr();

    /// Copy constructor.
    //  Transmits ownership.
    OwnPtr(const OwnPtr &);

    /// Constructor.
    //  Construct pointer to an object just created by "new". Grabs ownership.
    OwnPtr(Object *);

    /// Destructor.
    //  Delete the pointed-at object, if the pointer is non-null.
    ~OwnPtr();

    /// Assign.
    //  Transmit ownership to copy.
    OwnPtr &operator=(const OwnPtr &);

    /// Assign.
    //  Grab the ownership of an object just created by "new".
    OwnPtr &operator=(Object *);

    /// Delegation operator.
    //  Return a C-type pointer to the object.
    Object *operator->() const;

    /// Dereferencing operator.
    //  Return a reference to the object.
    Object &operator*() const;

    /// Test for validity.
    //  Return true, if the pointer is not NULL.
    bool isValid() const;

    /// Release ownership.
    //  Return built-in pointer.
    //  The calling program must take over the ownership.
    Object *release();

private:

    // The object pointed at.
    mutable Object *object;
};


// Inline functions for template class OwnPtr
// ------------------------------------------------------------------------

template <class Object>
inline OwnPtr<Object>::OwnPtr():
    object(0)
{}


template <class Object>
inline OwnPtr<Object>::OwnPtr(const OwnPtr<Object> &rhs):
    object(rhs.object) {
    rhs.object = 0;
}


template <class Object>
inline OwnPtr<Object>::OwnPtr(Object *obj):
    object(obj)
{}


template <class Object>
inline OwnPtr<Object>::~OwnPtr() {
    delete object;
}


template <class Object>
inline OwnPtr<Object> &OwnPtr<Object>::operator=(const OwnPtr &rhs) {
    if(object != rhs.object) {
        delete object;
        object = rhs.object;
        rhs.object = 0;
    }

    return *this;
}


template <class Object>
inline OwnPtr<Object> &OwnPtr<Object>::operator=(Object *obj) {
    delete object;
    object = obj;
    return *this;
}


template <class Object>
inline Object *OwnPtr<Object>::operator->() const {
    return object;
}


template <class Object>
inline Object &OwnPtr<Object>::operator*() const {
    return *object;
}


template <class Object>
inline bool OwnPtr<Object>::isValid() const {
    return object != 0;
}


template <class Object>
inline Object *OwnPtr<Object>::release() {
    Object *temp = object;
    object = 0;
    return temp;
}

#endif // CLASSIC_OwnPtr_HH
