// ------------------------------------------------------------------------
// $RCSfile: IfStatement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: IfStatement
//   Representation for OPAL IF statements.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/IfStatement.h"

#include "OpalParser/CompoundStatement.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Parser/Parser.h"
#include "Parser/Token.h"
#include "Parser/TokenStream.h"
#include "Utilities/ParseError.h"

#include "Utility/IpplInfo.h"

// class IfStatement
//   Statement of the form "IF ( <condition> ) <statement>".
// ------------------------------------------------------------------------

IfStatement::IfStatement(const Parser &parser, TokenStream &is):
    Statement("", 0), then_block(0), else_block(0) {
    Token key = is.readToken();
    Token token = is.readToken();

    if(key.isKey("IF") && token.isDel('(')) {
        append(token);
        token = is.readToken();
        int level = 1;

        while(! token.isEOF()) {
            append(token);

            if(token.isDel('(')) {
                level++;
            } else if(token.isDel(')')) {
                level--;
                if(level == 0) break;
            }

            token = is.readToken();
        }

        then_block = parser.readStatement(&is);
        token = is.readToken();

        if(! token.isEOF() && token.isKey("ELSE")) {
            else_block = parser.readStatement(&is);
        } else {
            is.putBack(token);
        }
    } else {
        throw ParseError("IfStatement::IfStatement()",
                         "Invalid \"IF\" statement.");
    }
}


IfStatement::~IfStatement() {
    delete then_block;
    delete else_block;
}


void IfStatement::execute(const Parser &parser) {
    start();
    Attribute condition = Attributes::makeBool("IF()", "");

    try {
        condition.parse(*this, false);
        OpalData::getInstance()->update();

        if(Attributes::getBool(condition)) {
            then_block->execute(parser);
        } else if(else_block) {
            else_block->execute(parser);
        }
    } catch(...) {
        std::ostringstream oss;
        this->print(oss);
        ERRORMSG("Invalid IF condition '" + oss.str() + "'");
        throw;
    }
}