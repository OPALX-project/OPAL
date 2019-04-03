#ifndef OPAL_Attribute_HH
#define OPAL_Attribute_HH

// ------------------------------------------------------------------------
// $RCSfile: Attribute.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Attribute:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#include "AbstractObjects/AttributeHandler.h"
#include "MemoryManagement/Pointer.h"
#include <iosfwd>
#include <string>

// #include <Ippl.h>
#include <Utility/Inform.h>

class Statement;

// Class Attribute
// ------------------------------------------------------------------------
/// A representation of an Object attribute.
//  Contains a pointer [b]base[/b] to the polymorphic attribute value,
//  derived from [b]AttributeBase[/b], and a pointer [b]handler[/b] to
//  the polymorphic attribute parser, derived from [b]AttributeHandler[/b].
//
//  The type of the attribute is determined by the class of AttributeHandler
//  attached.  The attribute handler must be set for all attributes by the
//  exemplar constructor for a command object.  The clone constructor does
//  not copy the handlers, but it creates new references to the same handlers.
//
//  The value is stored in the object pointed at by base.  An object can be
//  created for the default value of an attribute by the exemplar constructor,
//  or it can be null, meaning an as yet undefined value.  The validity of
//  the pointer can be tested by the method Attribute::isValid().
//
//  Each OPAL command object (derived from class Object) has a
//  vector<Attribute> containing the command attributes read by the parser.
//  This attribute vector is first filled by the clone constructor of the
//  command class by copying the pointers in the model, and possibly by the
//  OPAL parser with the values given by the user.

class Attribute {

public:

    /// Default constructor.
    //  Leaves both pointers [b]NULL[/b].
    //  An object constructed by the default constructor must be initialized
    //  using assignment before it can be used.
    //  Failing to do so may cause a program crash.
    Attribute();

    /// Copy constructor.
    //  Both value and parser are shared with [b]rhs[/b].
    Attribute(const Attribute &rhs);

    /// Constructor defining a parser and an optional value.
    //  The default value is ``undefined''.
    Attribute(AttributeHandler *h, AttributeBase *b = nullptr);

    ~Attribute();
    const Attribute &operator=(const Attribute &);

    /// Test for valid pointer.
    //  Returns true, when the value is defined (non-null pointer).
    operator bool() const;

    /// Return reference to polymorphic value.
    AttributeBase &getBase() const;

    bool isBaseAllocated() const;

    /// Return printable representation.
    //  This string duplicates the input expression defining the attribute.
    std::string getImage() const;

    /// Return a reference to the parser.
    AttributeHandler &getHandler() const;

    /// Return the help string.
    //  This string is stored in the parser object.
    const std::string &getHelp() const;

    /// Return the attribute name.
    //  This string is stored in the parser object.
    const std::string &getName() const;

    /// Return the attribute type.
    //  This string ("real", "logical", etc.) is stored in the parser object.
    const std::string &getType() const;

    /// Return [b]deferred[/b] flag.
    //  If this flag is set, any expression is re-evaluated each time the
    //  value is fetched.  Normally attribute evaluation is not deferred,
    //  i. e. any expression is evaluated only the first time the value is
    //  fetched after a new definition has been read.
    //  See [b]Expressions::ADeferred[/b] and [b]Expressions::SDeferred[/b]
    //  for information.
    bool isDeferred() const;

    /// Set read-only flag.
    //  If this flag is set, the attribute cannot be redefined by the user.
    //  See [b]Expressions::ADeferred[/b] and [b]Expressions::SDeferred[/b]
    //  for information.
    void setDeferred(bool);

    /// Test for expression.
    //  Return true, if the attribute is defined as an expression.
    bool isExpression() const;

    /// Test for read only.
    //  Return true, if the attribute cannot be redefined by the user.
    bool isReadOnly() const;

    /// Set read-only flag.
    //  Set or reset the attribute's read-only flag.
    void setReadOnly(bool);

    /// Parse attribute.
    //  Use the contained parser to parse a new value for the attribute.
    void parse(Statement &stat, bool eval);

    /// Parse array component.
    //  Use the contained parser to parse a new value for an existing
    //  vector component.  If the attribute is a scalar value, or if the
    //  component does not exist, this method throws [b]OpalException[/b].
    void parseComponent(Statement &stat, bool eval, int index);

    /// Define new value.
    //  Assign a new value.  The value must be compatible with the parser
    //  type, otherwise a [b]OpalException[/b] is thrown.
    void set(AttributeBase *newBase);

    /// Assign default value.
    //  Set the attribute value to the default value stored in the parser.
    //  If no default value exists, set it to ``undefined''.
    void setDefault();

    /// Print attribute.
    //  Print the attribute name, followed by an equals sign and its value.
    void print(int &len) const;

    bool defaultUsed() const;
private:

    // Pointer to the value.  The value can be shared for several objects.
    Pointer<AttributeBase> base;

    // Pointer to the shared attribute parser.
    Pointer<AttributeHandler> handler;

    bool isDefault;
};


// Class Attribute
// ------------------------------------------------------------------------

inline Attribute::operator bool() const {
    return base != nullptr;
}

inline bool Attribute::defaultUsed() const {
    return isDefault;
}

extern std::ostream &operator<<(std::ostream &os, const Attribute &attr);

#endif // OPAL_Attribute_HH