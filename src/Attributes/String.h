#ifndef OPAL_String_HH
#define OPAL_String_HH

// ------------------------------------------------------------------------
// $RCSfile: String.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: String
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"


// Class String
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type string.
    class String: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        String(const std::string &name, const std::string &help);

        virtual ~String();

        /// Return attribute type string ``string''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        String();
        String(const String &);
        void operator=(const String &);
    };

};

#endif // OPAL_String_HH
