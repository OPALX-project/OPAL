#ifndef CLASSIC_ConstChannel_HH
#define CLASSIC_ConstChannel_HH

// ------------------------------------------------------------------------
// $RCSfile: ConstChannel.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstChannel
//
// ------------------------------------------------------------------------
// Class category: Channels
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Class ConstChannel
// ------------------------------------------------------------------------
/// Abstract interface for read-only access to variable.
//  Class ConstChannel allows access to an arbitrary [b]double[/b] value.

class ConstChannel {

public:

    ConstChannel();
    virtual ~ConstChannel();

    /// Duplicate the channel.
    virtual ConstChannel *clone() const = 0;

    /// Read channel.
    //  If the channel can be read, set [b]value[/b] and return true,
    //  otherwise return false.
    virtual bool get(double &value) const = 0;

    /// Read channel.
    //  Return the value read from the channel.
    operator double() const
    { double value = 0.0; get(value); return value; }

    /// Check if settable.
    //  Return false for this class (channel cannot be set).
    virtual bool isSettable() const;

private:

    // Not implemented.
    ConstChannel(const ConstChannel &);
    const ConstChannel &operator=(const ConstChannel &);
};

#endif // CLASSIC_ConstChannel_HH
