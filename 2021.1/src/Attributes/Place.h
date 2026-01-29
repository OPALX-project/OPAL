#ifndef OPAL_Place_HH
#define OPAL_Place_HH

// ------------------------------------------------------------------------
// $RCSfile: Place.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Place
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/PlaceRep.h"

class Attribute;

// Class Place
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type place reference.
    class Place: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Place(const std::string &name, const std::string &help);

        virtual ~Place();

        /// Return attribute type string ``place''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        Place();
        Place(const Place &);
        void operator=(const Place &);
    };

};

#endif // OPAL_Place_HH
