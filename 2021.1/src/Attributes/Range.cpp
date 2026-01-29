// ------------------------------------------------------------------------
// $RCSfile: Range.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Range
//   A class used to parse range attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/Range.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using namespace Expressions;


// Class Range
// ------------------------------------------------------------------------

namespace Attributes {

    Range::Range(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    Range::~Range()
    {}


    const std::string &Range::getType() const {
        static const std::string type("range");
        return type;
    }


    void Range::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue<RangeRep>(parseRange(stat)));
    }

};
