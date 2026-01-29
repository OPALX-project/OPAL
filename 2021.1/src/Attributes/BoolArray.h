#ifndef OPAL_BoolArray_HH
#define OPAL_BoolArray_HH

// ------------------------------------------------------------------------
// $RCSfile: BoolArray.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class BoolArray:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"


// Class BoolArray
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type logical array.
    class BoolArray: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        BoolArray(const std::string &name, const std::string &help);

        virtual ~BoolArray();

        /// Return attribute type string ``logical array''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

        /// Parse a component of the array.
        //  Identified by its index.
        virtual void parseComponent(Attribute &, Statement &, bool, int) const;

    private:

        // Not implemented.
        BoolArray();
        BoolArray(const BoolArray &);
        void operator=(const BoolArray &);
    };

};

#endif // OPAL_BoolArray_HH
