// ------------------------------------------------------------------------
// $RCSfile: TokenList.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class TokenList:
//   A class used to parse a "token list" attribute.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/TokenList.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using namespace Expressions;


// Class TokenList
// ------------------------------------------------------------------------

namespace Attributes {

    TokenList::TokenList(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    TokenList::~TokenList()
    {}


    const std::string &TokenList::getType() const {
        static const std::string type("token list");
        return type;
    }


    void TokenList::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue<std::list<Token> >(parseTokenList(stat)));
    }

};
