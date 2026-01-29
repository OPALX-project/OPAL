#ifndef OPAL_WhileStatement_HH
#define OPAL_WhileStatement_HH 1

// ------------------------------------------------------------------------
// $RCSfile: WhileStatement.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: WhileStatement
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Parser/Statement.h"
#include <iosfwd>

class Parser;
class TokenStream;


// class "WhileStatement":
// ------------------------------------------------------------------------
/// While statement.
//  A statement of the form "WHILE ( <condition> ) <statement>".
//  The condition is stored in the Token list inherited from Statement,
//  the block to be executed repeatedly in a Statement.

class WhileStatement: public Statement {

public:

    /// Constructor.
    //  Parse the statement on the given token stream, using the given parser.
    WhileStatement(const Parser &, TokenStream &);

    virtual ~WhileStatement();

    /// Execute.
    //  Use the given parser to execute the controlled statements.
    virtual void execute(const Parser &);

private:

    // Not implemented.
    WhileStatement();
    WhileStatement(const WhileStatement &);
    void operator=(const WhileStatement &);

    Statement *while_block;
};

#endif // OPAL_WhileStatement_HH
