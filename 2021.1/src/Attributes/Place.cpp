// ------------------------------------------------------------------------
// $RCSfile: Place.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Place
//   A class used to parse place attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/Place.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using namespace Expressions;


// Class Place
// ------------------------------------------------------------------------

namespace Attributes {

    Place::Place(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    Place::~Place()
    {}


    const std::string &Place::getType() const {
        static const std::string type("place");
        return type;
    }


    void Place::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue<PlaceRep>(parsePlace(stat)));
    }

};
