#ifndef OPAL_IfStatement_HH
#define OPAL_IfStatement_HH 1

// ------------------------------------------------------------------------
// $RCSfile: IfStatement.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: IfStatement
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Parser/Statement.h"
#include "AbstractObjects/Attribute.h"

class Parser;
class TokenStream;


// class "IfStatement":
// ------------------------------------------------------------------------
/// If statement.
//  A statement of the form "IF ( <condition> ) <statement>".
//  The condition is stored in the Token list inherited from Statement,
//  the then and else blocks are stored in two Statements.

class IfStatement: public Statement {

public:

    /// Constructor.
    //  Parse the statement on the given token stream, using the given parser.
    IfStatement(const Parser &, TokenStream &);

    virtual ~IfStatement();

    /// Execute.
    //  Use the given parser to execute the controlled statements.
    virtual void execute(const Parser &);

private:

    // Not implemented.
    IfStatement();
    IfStatement(const IfStatement &);
    void operator=(const IfStatement &);

    Statement *then_block;
    Statement *else_block;
};

#endif // OPAL_IfStatement_HH
