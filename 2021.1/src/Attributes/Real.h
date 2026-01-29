#ifndef OPAL_Real_HH
#define OPAL_Real_HH

// ------------------------------------------------------------------------
// $RCSfile: Real.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Real
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"


// Class Real
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type real.
    class Real: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Real(const std::string &name, const std::string &help);

        virtual ~Real();

        /// Return attribute type string ``real''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        Real();
        Real(const Real &);
        void operator=(const Real &);
    };

};

#endif // OPAL_Real_HH
