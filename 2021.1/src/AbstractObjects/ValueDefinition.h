#ifndef OPAL_ValueDefinition_HH
#define OPAL_ValueDefinition_HH

// ------------------------------------------------------------------------
// $RCSfile: ValueDefinition.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ValueDefinition
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Object.h"


// Class ValueDefinition
// ------------------------------------------------------------------------
/// The base class for all OPAL value definitions.
//  It implements the common behaviour of scalar and vector variables,
//  it can also be used via dynamic casting to determine whether an object
//  represents a variable.

class ValueDefinition: public Object {

public:

    virtual ~ValueDefinition();

    /// Return the object category as a string.
    //  Return the string "VARIABLE".
    virtual const std::string getCategory() const;

    /// Trace flag.
    //  If true, the object's execute() function should be traced.
    //  Always false for value definitions.
    virtual bool shouldTrace() const;

    /// Update flag.
    //  If true, the data structure should be updated before calling execute().
    //  Always false for value definitions.
    virtual bool shouldUpdate() const;


    /// Return logical value.
    //  The default version throws OpalException.
    virtual bool getBool() const;

    /// Return indexed logical value.
    //  The default version throws OpalException.
    virtual bool getBoolComponent(int) const;

    /// Return real value.
    //  The default version throws OpalException.
    virtual double getReal() const;

    /// Return indexed real value.
    //  The default version throws OpalException.
    virtual double getRealComponent(int) const;

    /// Return string value.
    //  The default version throws OpalException.
    virtual std::string getString() const;

    /// Return indexed string value.
    //  The default version throws OpalException.
    virtual std::string getStringComponent(int) const;

    /// Return the attribute representing the value of the definition.
    //  Version for non-constant object.
    Attribute &value();

    /// Return the attribute representing the value of the definition.
    //  Version for constant object.
    const Attribute &value() const;

protected:

    /// Constructor for exemplars.
    ValueDefinition(int size, const char *name, const char *help);

    /// Constructor for clones.
    ValueDefinition(const std::string &name, ValueDefinition *parent);

    /// The declaration prefix.
    //  A string representing the value type in OPAL-9 input format:
    //  "BOOL", "REAL", "REAL CONST", "REAL VECTOR", etc.
    const std::string itsPrefix;

private:

    // Not implemented.
    ValueDefinition();
    ValueDefinition(const ValueDefinition &);
    void operator=(const ValueDefinition &);
};

#endif // OPAL_ValueDefinition_HH
