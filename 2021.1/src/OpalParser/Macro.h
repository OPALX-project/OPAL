#ifndef CLASSIC_Macro_HH
#define CLASSIC_Macro_HH

// ------------------------------------------------------------------------
// $RCSfile: Macro.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Macro
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Object.h"
#include "Parser/Token.h"
#include <vector>

class Statement;
class TokenStream;


// Class Macro
// ------------------------------------------------------------------------
/// Abstract base class for macros.
//  The base class for storing the ``archetype'' for all macro-like commands.

class Macro: public Object {

public:

    /// Exemplar constructor.
    Macro(int size, const char *name, const char *help);

    /// Clone constructor.
    Macro(const std::string &name, Object *parent);

    virtual ~Macro();

    /// Make clone.
    //  Throw ParseError, since for macro we must make a template and not
    //  a clone.
    virtual Macro *clone(const std::string &name);

    /// Return the object category as a string.
    //  Return the string "MACRO".
    virtual const std::string getCategory() const;

    /// Trace flag.
    //  If true, the object's execute() function should be traced.
    //  Always false for macros.
    virtual bool shouldTrace() const;

    /// Update flag.
    //  If true, the data structure should be updated before calling execute().
    //  Always false for macros.
    virtual bool shouldUpdate() const;


    /// Parse actual arguments.
    //  Expects parse pointer in the statement to be set on the first argument.
    virtual void parseActuals(Statement &);

    /// Parse formal arguments.
    //  Expects parse pointer in the statement to be set on the first argument.
    virtual void parseFormals(Statement &);

protected:

    /// The formal argument list.
    //  Each formal argument is represented as a name and stored as a string.
    std::vector<std::string> formals;

    /// The actual argument list.
    //  Each actual argument is stored as a vector of tokens.
    std::vector<std::vector < Token > > actuals;

private:

    // Not implemented.
    Macro();
    Macro(const Macro &);
    void operator=(const Macro &);
};

#endif // CLASSIC_Macro_HH
