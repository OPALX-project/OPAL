#ifndef OPAL_Range_HH
#define OPAL_Range_HH

// ------------------------------------------------------------------------
// $RCSfile: Range.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Range
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/RangeRep.h"

class Attribute;

// Class Range
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type range definition.
    class Range: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Range(const std::string &name, const std::string &help);

        virtual ~Range();

        /// Return attribute type ``range''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        Range();
        Range(const Range &);
        void operator=(const Range &);
    };

};

#endif // OPAL_Range_HH
