#ifndef CLASSIC_Parser_HH
#define CLASSIC_Parser_HH

// ------------------------------------------------------------------------
// $RCSfile: Parser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Parser
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

class Statement;
class TokenStream;


// Class Parser
// ------------------------------------------------------------------------
/// Interface for abstract language parser.

class Parser {

public:

    Parser();
    virtual ~Parser();

    /// Parse and execute the statement.
    virtual void parse(Statement &stat) const = 0;

    /// Read complete statement from token stream
    virtual Statement *readStatement(TokenStream *ts) const = 0;

    /// Read statements and parse.
    //  Read one statement at a time on token stream, parse it, and execute it.
    virtual void run(TokenStream *ts) const = 0;

private:

    // Not implemented.
    Parser(const Parser &);
    void operator=(const Parser &);
};

inline
Parser::Parser()
{ }

inline
Parser::~Parser()
{ }

#endif // CLASSIC_Parser_HH