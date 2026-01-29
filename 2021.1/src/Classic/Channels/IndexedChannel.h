#ifndef CLASSIC_IndexedChannel_HH
#define CLASSIC_IndexedChannel_HH

// ------------------------------------------------------------------------
// $RCSfile: IndexedChannel.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: IndexedChannel
//
// ------------------------------------------------------------------------
// Class category: Channels
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Channels/Channel.h"


// Class IndexedChannel
// ------------------------------------------------------------------------
/// Access to an indexed [b]double[/b] data member.
//  Template class IndirectChannel allows access to an indexed [b]double[/b]
//  data member of an object.

template <class T> class IndexedChannel: public Channel {

public:

    /// Constructor.
    //  The constructed channel provides access to an array of an object
    //  of class [b]T[/b].  The channel keeps a reference to [b]object[/b],
    //  the pointers to member [b]getF[/b] and [b]setF[/b], and the index
    //  [b]index[/b].
    //  [p]
    //  Values set are transmitted via [tt] object.*setF(index,value)[/tt]
    //  and read via [tt]value = object.*getF(index)[/tt].
    IndexedChannel(T &object, double(T::*getF)(int) const,
                   void (T::*setF)(int, double), int index);

    IndexedChannel(const IndexedChannel &);
    virtual ~IndexedChannel();

    /// Duplicate the channel.
    virtual IndexedChannel *clone() const;

    /// Fetch from channel.
    //  If the channel can be read, set [b]value[/b] and return true,
    //  otherwise return false.
    virtual bool get(double &) const;

    /// Store into channel.
    //  If the channel can be written, store [b]value[/b] into it and return true,
    //  otherwise return false.
    virtual bool set(double);

    /// Test if settable.
    //  Return true, if the channel can receive values, i.e. if the [b]sefF[/b]
    //  pointer is not NULL.
    virtual bool isSettable() const;

private:

    // Not implemented.
    IndexedChannel();
    const IndexedChannel &operator=(const IndexedChannel &);

    // Reference to the object to be set.
    T &itsObject;

    // The get and set functions for the channel.
    double(T::*getF)(int) const;
    void (T::*setF)(int, double);

    // The bias to be transmitted to the set() or get() method.
    int bias;
};


template <class T>
IndexedChannel<T>::IndexedChannel(T &object, double(T::*get)(int) const,
                                  void (T::*set)(int, double), int index):
    itsObject(object), getF(get), setF(set), bias(index)
{}


template <class T>
IndexedChannel<T>::IndexedChannel(const IndexedChannel &rhs):
    Channel(),
    itsObject(rhs.itsObject), getF(rhs.getF), setF(rhs.setF), bias(rhs.bias)
{}


template <class T>
IndexedChannel<T>::~IndexedChannel()
{}


template <class T>
IndexedChannel<T> *IndexedChannel<T>::clone() const {
    return new IndexedChannel(*this);
}


template <class T>
bool IndexedChannel<T>::get(double &value) const {
    value = (itsObject.*getF)(bias);
    return true;
}


template <class T>
bool IndexedChannel<T>::set(double value) {
    if(setF != 0) {
        (itsObject.*setF)(bias, value);
        return true;
    } else {
        return false;
    }
}


template <class T>
bool IndexedChannel<T>::isSettable() const {
    return (setF != 0);
}

#endif // CLASSIC_IndexedChannel_HH
