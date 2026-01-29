// ------------------------------------------------------------------------
// $RCSfile: SimpleStatement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SimpleStatement
//   Concrete representation for a standard input language statement.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Parser/SimpleStatement.h"
#include "Parser/Parser.h"
#include "Parser/Token.h"
#include "Parser/TokenStream.h"


// class SimpleStatement
// ------------------------------------------------------------------------

SimpleStatement::SimpleStatement(const std::string &name, int line):
    Statement(name, line)
{}


SimpleStatement::SimpleStatement(const std::string &name, TokenList &list):
    Statement(name, list)
{}


SimpleStatement::~SimpleStatement()
{}


void SimpleStatement::execute(const Parser &parser) {
    curr = keep = tokens.begin();
    if(curr != tokens.end()) parser.parse(*this);
}
