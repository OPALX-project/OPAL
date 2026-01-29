#ifndef OPAL_Bool_HH
#define OPAL_Bool_HH

// ------------------------------------------------------------------------
// $RCSfile: Bool.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Attributes::Bool
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"


// Class Attributes::Bool
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for attribute of type logical.
    class Bool: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Bool(const std::string &name, const std::string &help);

        virtual ~Bool();

        /// Return attribute type string ``logical''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        Bool();
        Bool(const Bool &);
        void operator=(const Bool &);
    };

};

#endif // OPAL_Bool_HH
