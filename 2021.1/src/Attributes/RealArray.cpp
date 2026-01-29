// ------------------------------------------------------------------------
// $RCSfile: RealArray.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class RealArray:
//   A class used to parse real array attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/RealArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AAutomatic.h"
#include "Expressions/ADeferred.h"
#include "Expressions/AList.h"
#include "Expressions/SConstant.h"
#include "Utilities/OpalException.h"
#include <vector>

using namespace Expressions;


// Class RealArray
// ------------------------------------------------------------------------

namespace Attributes {

    RealArray::RealArray(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    RealArray::~RealArray()
    {}


    const std::string &RealArray::getType() const {
        static std::string type = "real array";
        return type;
    }


    void RealArray::parse(Attribute &attr, Statement &statement, bool eval) const {
        PtrToArray<double> expr = parseRealArray(statement);

        if(eval) {
            // Use ADeferred here, since a component may be overridden later
            // by an expression.
            attr.set(new ADeferred<double>(expr->evaluate()));
        } else if(is_deferred) {
            attr.set(new ADeferred<double>(expr));
        } else {
            attr.set(new AAutomatic<double>(expr));
        }
    }


    void RealArray::parseComponent
    (Attribute &attr, Statement &stat, bool eval, int index) const {
        ADeferred<double> *array = 0;

        if(AttributeBase *base = &(attr.getBase())) {
            array = dynamic_cast<ADeferred<double>*>(base);
        } else {
            attr.set(array = new ADeferred<double>());
        }

        PtrToScalar<double> expr = parseReal(stat);
        if(eval) expr = new SConstant<double>(expr->evaluate());
        array->setComponent(index, expr);
    }

};
