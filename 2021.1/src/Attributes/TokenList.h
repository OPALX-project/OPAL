#ifndef OPAL_TokenList_HH
#define OPAL_TokenList_HH

// ------------------------------------------------------------------------
// $RCSfile: TokenList.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TokenList
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

template <class T> class TValue;


// Class TokenList
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type token list.
    //  Such an attribute can encode e.g. an expression in a LIST command.
    class TokenList: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        TokenList(const std::string &name, const std::string &help);

        virtual ~TokenList();

        /// Return attribute type string ``token list''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        TokenList();
        TokenList(const TokenList &);
        void operator=(const TokenList &);
    };

};

#endif // OPAL_TokenList_HH
