// ------------------------------------------------------------------------
// $RCSfile: Real.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Real:
//   A class used to parse real attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/12/15 10:08:28 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Attributes/Real.h"
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


// Class Real
// ------------------------------------------------------------------------

namespace Attributes {

    Real::Real(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    Real::~Real()
    {}


    const std::string &Real::getType() const {
        static const std::string type("real");
        return type;
    }


    void Real::parse(Attribute &attr, Statement &statement, bool eval) const {
        PtrToScalar<double> expr = parseReal(statement);

        if(eval || expr->isConstant()) {
            attr.set(new SValue<double>(expr->evaluate()));
        } else if(isDeferred()) {
            attr.set(new SDeferred<double>(expr));
        } else {
            attr.set(new SAutomatic<double>(expr));
        }
    }

};
