#ifndef OPAL_StringArray_HH
#define OPAL_StringArray_HH

// ------------------------------------------------------------------------
// $RCSfile: StringArray.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class StringArray:
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


// Class StringArray
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type string array.
    class StringArray: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        StringArray(const std::string &name, const std::string &help);

        virtual ~StringArray();

        /// Return attribute type string ``string array''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

        /// Parse a component of the array.
        //  Identified by its index.
        virtual void parseComponent(Attribute &, Statement &, bool, int) const;

    private:

        // Not implemented.
        StringArray();
        StringArray(const StringArray &);
        void operator=(const StringArray &);
    };

};

#endif // OPAL_StringArray_HH
