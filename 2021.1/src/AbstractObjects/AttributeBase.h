#ifndef OPAL_AttributeBase_HH
#define OPAL_AttributeBase_HH

// ------------------------------------------------------------------------
// $RCSfile: AttributeBase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeBase
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:40:33 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "MemoryManagement/RCObject.h"
#include "AbstractObjects/Invalidator.h"
#include <iosfwd>
#include <string>

// Class AttributeBase
// ------------------------------------------------------------------------
/// Abstract base class for attribute values of different types.
//  The various derived classes contain the actual values, notably the
//  template classes SValue<T> for scalars and AValue<T> for arrays.

class AttributeBase: public RCObject, public Invalidator {

public:

    AttributeBase();
    virtual ~AttributeBase();

    /// Make clone.
    //  Construct an exact copy of the value.
    virtual AttributeBase *clone() const = 0;

    /// Convert to string.
    //  Uses [b]print()[/b] to convert the input expression for the attribute
    //  to a printable representation.
    std::string getImage() const;

    /// Test for expression.
    //  Return [b]true[/b] if the value is an expression.
    virtual bool isExpression() const;

    /// Print value.
    //  Print the value on the given output stream.
    //  The result allows reconstruction of any expression.
    virtual void print(std::ostream &) const = 0;

private:

    // Not implemented.
    AttributeBase(AttributeBase &);
    void operator=(AttributeBase &);
};

#endif // OPAL_AttributeBase_HH
