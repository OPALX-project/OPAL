// ------------------------------------------------------------------------
// $RCSfile: Bool.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Attributes::Bool
//   A class used to parse boolean attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/Bool.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SDeferred.h"
#include "Expressions/SValue.h"
#include "Parser/SimpleStatement.h"
#include "Parser/StringStream.h"
#include "Parser/Token.h"
#include "Utilities/ClassicException.h"
#include "Utilities/ParseError.h"

using namespace Expressions;
using std::cerr;
using std::endl;


// Class Attributes::Bool
// ------------------------------------------------------------------------

namespace Attributes {

    Bool::Bool(const std::string &name, const std::string &help):
        AttributeHandler(name, help, new SValue<bool>(true))
    {}


    Bool::~Bool()
    {}

    const std::string &Bool::getType() const {
        static const std::string type("logical");
        return type;
    }

    void Bool::parse(Attribute &attr, Statement &stat, bool eval) const {
        PtrToScalar<bool> expr = parseBool(stat);

        if(eval || expr->isConstant()) {
            attr.set(new SValue<bool>(expr->evaluate()));
        } else if(is_deferred) {
            attr.set(new SDeferred<bool>(expr));
        } else {
            attr.set(new SAutomatic<bool>(expr));
        }
    }

};
