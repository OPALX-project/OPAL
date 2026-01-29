#ifndef OPAL_TokenListArray_HH
#define OPAL_TokenListArray_HH

// ------------------------------------------------------------------------
// $RCSfile: TokenListArray.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TokenListArray
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"
#include "Parser/Token.h"
#include <list>


// Class TokenListArray
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type token list array.
    //  Such an attribute may encode a list of expressions in a LIST command.
    class TokenListArray: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        TokenListArray(const std::string &name, const std::string &help);

        virtual ~TokenListArray();

        /// Return attribute type string ``token list array''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

        /// Parse a component.
        //  Identified by its index.
        virtual void parseComponent(Attribute &, Statement &, bool, int) const;

    private:

        // Not implemented.
        TokenListArray();
        TokenListArray(const TokenListArray &);
        void operator=(const TokenListArray &);
    };

};

#endif // OPAL_TokenListArray_HH
