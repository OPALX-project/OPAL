#ifndef CLASSIC_SimpleStatement_HH
#define CLASSIC_SimpleStatement_HH 1

// ------------------------------------------------------------------------
// $RCSfile: SimpleStatement.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Statement
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Parser/Statement.h"


// class SimpleStatement
// ------------------------------------------------------------------------
/// A simple input statement in token form.
//  The statement is stored as a list of Token's.

class SimpleStatement: public Statement {

public:

    /// Constructor.
    //  Stores the name of the input stream and the line at which
    //  the statement begins.
    SimpleStatement(const std::string &streamName, int streamLine);

    /// Constructor.
    //  Stores a name (e.g. for a macro) and the token list.
    SimpleStatement(const std::string &streamName, TokenList &list);

    /// Destructor.
    virtual ~SimpleStatement();

    /// Execute statement.
    //  This method provides a hook for structured statement execution.
    //  It is called by the parser, and may be overridden to execute
    //  e.g. a conditional or a loop statement.
    virtual void execute(const Parser &);
};

#endif // CLASSIC_SimpleStatement_HH
