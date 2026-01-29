// ------------------------------------------------------------------------
// $RCSfile: Reference.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Reference:
//   A class used to parse reference attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/Reference.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SRefAttr.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using namespace Expressions;


// Class Reference
// ------------------------------------------------------------------------

namespace Attributes {

    Reference::Reference(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    Reference::~Reference()
    {}


    const std::string &Reference::getType() const {
        static const std::string type("reference");
        return type;
    }


    void Reference::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue < SRefAttr<double> > (*parseReference(stat)));
    }

};
