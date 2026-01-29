#ifndef OPAL_Reference_HH
#define OPAL_Reference_HH

// ------------------------------------------------------------------------
// $RCSfile: Reference.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Reference:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"

class Attribute;


// Class Reference
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type attribute reference.
    //  The attribute referred to may be logical, real, or string.
    class Reference: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Reference(const std::string &name, const std::string &help);

        virtual ~Reference();

        /// Return attribute type string ``reference''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        Reference(const Reference &);
        void operator=(const Reference &);
    };

};

#endif // OPAL_Reference_HH
