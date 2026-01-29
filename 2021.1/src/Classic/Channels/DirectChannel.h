#ifndef CLASSIC_DirectChannel_HH
#define CLASSIC_DirectChannel_HH

// ------------------------------------------------------------------------
// $RCSfile: DirectChannel.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DirectChannel
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


// Class DirectChannel
// ------------------------------------------------------------------------
/// Direct access to a [b]double[/b] variable.
//  Class DirectChannel allows direct access to a [b]double[/b] variable.

class DirectChannel: public Channel {

public:

    /// Constructor.
    //  The constructed channel gives read/write access to the variable
    //  [b]value[/b].  The variable [b]value[/b] must not be destroyed as
    //  long as this channel is active.
    explicit DirectChannel(double &value);

    DirectChannel(const DirectChannel &);
    virtual ~DirectChannel();

    /// Duplicate the channel.
    virtual DirectChannel *clone() const;

    /// Store into channel.
    //  If the channel can be written,
    //  store [b]value[/b] into it and return true,
    //  otherwise return false.
    virtual bool set(double);

    /// Fetch from channel.
    //  If the channel can be read, set [b]value[/b] and return true,
    //  otherwise return false.
    virtual bool get(double &) const;

private:

    // Not implemented.
    DirectChannel();
    const DirectChannel &operator=(const DirectChannel &);

    // The address of the variable to be read or written.
    double &reference;
};


inline DirectChannel::DirectChannel(double &value):
    reference(value)
{}


inline DirectChannel::DirectChannel(const DirectChannel &rhs):
    Channel(),
    reference(rhs.reference)
{}


inline DirectChannel::~DirectChannel()
{}


inline DirectChannel *DirectChannel::clone() const {
    return new DirectChannel(*this);
}


inline bool DirectChannel::set(double value) {
    reference = value;
    return true;
}


inline bool DirectChannel::get(double &value) const {
    value = reference;
    return true;
}


#endif // CLASSIC_DirectChannel_HH
