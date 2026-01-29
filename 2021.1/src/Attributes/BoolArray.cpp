// ------------------------------------------------------------------------
// $RCSfile: BoolArray.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class BoolArray:
//   A class used to parse boolean array attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/BoolArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AAutomatic.h"
#include "Expressions/ADeferred.h"
#include "Expressions/AList.h"
#include "Expressions/SConstant.h"
#include "Utilities/OpalException.h"
#include <vector>

using namespace Expressions;


// Class BoolArray
// ------------------------------------------------------------------------

namespace Attributes {

    BoolArray::BoolArray(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    BoolArray::~BoolArray()
    {}


    const std::string &BoolArray::getType() const {
        static std::string type = "logical array";
        return type;
    }


    void BoolArray::parse(Attribute &attr, Statement &stat, bool eval) const {
        Expressions::PtrToArray<bool> expr = Expressions::parseBoolArray(stat);

        if(eval) {
            // Use ADeferred here, since a component may be overridden later
            // by an expression.
            attr.set(new ADeferred<bool>(expr->evaluate()));
        } else if(is_deferred) {
            attr.set(new ADeferred<bool>(expr));
        } else {
            attr.set(new AAutomatic<bool>(expr));
        }
    }


    void BoolArray::parseComponent
    (Attribute &attr, Statement &stat, bool eval, int index) const {
        ADeferred<bool> *array = 0;

        if(AttributeBase *base = &(attr.getBase())) {
            array = dynamic_cast<ADeferred<bool>*>(base);
        } else {
            attr.set(array = new ADeferred<bool>());
        }

        PtrToScalar<bool> expr = parseBool(stat);
        if(eval) expr = new SConstant<bool>(expr->evaluate());
        array->setComponent(index, expr);
    }

};
