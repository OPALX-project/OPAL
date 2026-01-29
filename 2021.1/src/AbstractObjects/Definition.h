#ifndef OPAL_Definition_HH
#define OPAL_Definition_HH

// ------------------------------------------------------------------------
// $RCSfile: Definition.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Definition
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Object.h"


// Class Definition
// ------------------------------------------------------------------------
/// The base class for all OPAL definitions.
//  It implements the common behaviour of definitions, it can also be used
//  via dynamic casting to determine whether an object represents a definition.

class Definition: public Object {

public:

    virtual ~Definition();

    /// Return the object category as a string.
    //  Return the string "DEFINITION".
    virtual const std::string getCategory() const;

    /// Trace flag.
    //  If true, the object's execute() function should be traced.
    //  Always false for definitions.
    virtual bool shouldTrace() const;

    /// Update flag.
    //  If true, the data structure should be updated before calling execute().
    //  Always false for definitions.
    virtual bool shouldUpdate() const;

protected:

    /// Constructor for exemplars.
    Definition(int size, const char *name, const char *help);

    /// Constructor for clones.
    Definition(const std::string &name, Definition *parent);

private:

    // Not implemented.
    Definition();
    Definition(const Definition &);
    void operator=(const Definition &);
};

#endif // OPAL_Definition_HH
