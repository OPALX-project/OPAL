#ifndef CLASSIC_AttributeSet_HH
#define CLASSIC_AttributeSet_HH

// ------------------------------------------------------------------------
// $RCSfile: AttributeSet.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeSet
//   A map of name (std::string) versus value (double) intended to store
//   user-defined attributes.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/12/16 16:26:43 $
// $Author: mad $
//
// ------------------------------------------------------------------------

//#include "Channels/Channel.h"
#include <map>
#include <functional>
#include <string>

class Channel;
class ConstChannel;

// Class AttributeSet
// ------------------------------------------------------------------------
/// Map of std::string versus double value.
//  Class AttributeSet implements a map of name (std::string) versus value
//  (double) for user-defined attributes.  This map is intended for
//  algorithms that require specific, but not predefined data in the
//  accelerator model for their working.

class AttributeSet {

public:

    /// A map of name versus value.
    typedef std::map<std::string, double, std::less<std::string> > NameMap;

    /// An iterator for a map of name versus value.
    typedef NameMap::const_iterator const_iterator;

    /// Default constructor.
    //  Constructs an empty map.
    AttributeSet();

    AttributeSet(const AttributeSet &);
    virtual ~AttributeSet();
    const AttributeSet &operator=(const AttributeSet &);

    /// Iterator accessing first member.
    const_iterator begin() const;

    /// Iterator marking the end of the list.
    const_iterator end() const;


    /// Get attribute value.
    //  If the attribute does not exist, return zero.
    double getAttribute(const std::string &aKey) const;

    /// Test for presence of an attribute.
    //  If the attribute exists, return true, otherwise false.
    bool hasAttribute(const std::string &aKey) const;

    /// Remove an existing attribute.
    //  If the key [b]aKey[/b] exists, this method removes it.
    void removeAttribute(const std::string &aKey);

    /// Set value of an attribute.
    void setAttribute(const std::string &aKey, double val);


    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    Channel *getChannel(const std::string &aKey, bool create = false);

    /// Construct a read-only channel.
    //  This method constructs a Channel permitting read-only access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    const ConstChannel *getConstChannel(const std::string &aKey) const;

protected:

    /// The attribute map.

    NameMap itsMap;
};


// Implementation.
// ------------------------------------------------------------------------

inline AttributeSet::const_iterator AttributeSet::begin() const
{ return itsMap.begin(); }

inline AttributeSet::const_iterator AttributeSet::end() const
{ return itsMap.end(); }

#endif // CLASSIC_AttributeSet_HH
