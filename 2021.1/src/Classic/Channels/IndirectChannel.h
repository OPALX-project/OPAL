#ifndef CLASSIC_IndirectChannel_HH
#define CLASSIC_IndirectChannel_HH

// ------------------------------------------------------------------------
// $RCSfile: IndirectChannel.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: IndirectChannel
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


// Class IndirectChannel
// ------------------------------------------------------------------------
/// Access to a [b]double[/b] data member.
//  Template class IndirectChannel allows access to a [b]double[/b] data
//  member of an object.

template <class T> class IndirectChannel: public Channel {

public:

    /// Constructor.
    //  The constructed channel provides access to a member of an object
    //  of class [b]T[/b].  The channel keeps a reference to [b]object[/b]
    //  and the pointers to member [b]getF[/b] and [b]setF[/b].
    //  Values set are transmitted via  object.*setF(value)
    //  and read via  value = object.*getF().
    IndirectChannel(T &object, double(T::*getF)() const,
                    void (T::*setF)(double));

    IndirectChannel(const IndirectChannel &);
    virtual ~IndirectChannel();

    /// Duplicate the channel.
    virtual IndirectChannel *clone() const;

    /// Fetch from channel.
    //  If the channel can be read, set [b]value[/b] and return true,
    //  otherwise return false.
    virtual bool get(double &) const;

    /// Store into channel.
    //  If the channel can be written,
    //  store [b]value[/b] into it and return true,
    //  otherwise return false.
    virtual bool set(double);

    /// Test if settable.
    //  Return true, if the channel can be written, i.e. if the set
    //  method pointer is not NULL
    virtual bool isSettable() const;

private:

    // Not implemented.
    IndirectChannel();
    const IndirectChannel &operator=(const IndirectChannel &);

    // Reference to the object to be set.
    T &itsObject;

    // The get and set functions for the channel.
    double(T::*getF)() const;
    void (T::*setF)(double);
};


template <class T>
IndirectChannel<T>::IndirectChannel(T &object, double(T::*get)() const,
                                    void (T::*set)(double)):
    itsObject(object), getF(get), setF(set)
{}


template <class T>
IndirectChannel<T>::IndirectChannel(const IndirectChannel &rhs):
    Channel(),
    itsObject(rhs.itsObject), getF(rhs.getF), setF(rhs.setF)
{}


template <class T>
IndirectChannel<T>::~IndirectChannel()
{}


template <class T>
IndirectChannel<T> *IndirectChannel<T>::clone() const {
    return new IndirectChannel(*this);
}


template <class T>
bool IndirectChannel<T>::get(double &value) const {
    value = (itsObject.*getF)();
    return true;
}


template <class T>
bool IndirectChannel<T>::set(double value) {
    if(setF != 0) {
        (itsObject.*setF)(value);
        return true;
    } else {
        return false;
    }
}


template <class T> inline
bool IndirectChannel<T>::isSettable() const {
    return (setF != 0);
}

#endif // CLASSIC_IndirectChannel_HH
