#ifndef OPAL_RealArray_HH
#define OPAL_RealArray_HH

// ------------------------------------------------------------------------
// $RCSfile: RealArray.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class RealArray:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"


// Class RealArray
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type real array.
    class RealArray: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        RealArray(const std::string &name, const std::string &help);

        virtual ~RealArray();

        /// Return attribute type string ``real array''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

        /// Parse a component of the array.
        //  Identified by its index.
        virtual void parseComponent(Attribute &, Statement &, bool, int) const;

    private:

        // Not implemented.
        RealArray();
        RealArray(const RealArray &);
        void operator=(const RealArray &);
    };

};

#endif // OPAL_RealArray_HH
