#ifndef CLASSIC_Channel_HH
#define CLASSIC_Channel_HH

// ------------------------------------------------------------------------
// $RCSfile: Channel.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Channel
//
// ------------------------------------------------------------------------
// Class category: Channels
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Channels/ConstChannel.h"


// Class Channel
// ------------------------------------------------------------------------
/// Abstract interface for read/write access to variable.
//  Class Channel allows access to an arbitrary [b]double[/b] value.


class Channel: public ConstChannel {

public:

    Channel();
    virtual ~Channel();

    /// Duplicate the channel.
    virtual Channel *clone() const = 0;

    /// Assign [b]value[/b] to channel.
    double operator=(double value)
    { set(value); return value; }

    /// Add and assign [b]value[/b] to channel.
    double operator+=(double value)
    { double temp; get(temp); temp += value; set(temp); return temp; }

    /// Subtract and assign [b]value[/b] to channel.
    double operator-=(double value)
    { double temp; get(temp); temp -= value; set(temp); return temp; }

    /// Store [b]value[/b] into channel.
    //  If the channel can be written,
    //  store [b]value[/b] into it and return true,
    //  otherwise return false.
    virtual bool set(double value) = 0;

    /// Test if settable.
    //  Return true, if the channel can receive values.
    virtual bool isSettable() const;

private:

    // Not implemented.
    Channel(const Channel &);
    const Channel &operator=(const Channel &);
};

#endif // CLASSIC_Channel_HH
