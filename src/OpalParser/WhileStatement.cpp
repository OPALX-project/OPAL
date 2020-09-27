// ------------------------------------------------------------------------
// $RCSfile: WhileStatement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: WhileStatement
//   Representation for OPAL WHILE statement.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/WhileStatement.h"

#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "OpalParser/CompoundStatement.h"
#include "Parser/Parser.h"
#include "Parser/Token.h"
#include "Parser/TokenStream.h"
#include "Utilities/ParseError.h"

#include "Utility/IpplInfo.h"

// class WhileStatement
//   Statement of the form "WHILE ( <condition> ) <statement>".
// ------------------------------------------------------------------------

WhileStatement::WhileStatement(const Parser &parser, TokenStream &is):
    Statement("", 0), while_block(0) {
    Token key = is.readToken();
    Token token = is.readToken();

    if(key.isKey("WHILE") && token.isDel('(')) {
        int level = 1;
        append(token);
        token = is.readToken();

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

        while_block = parser.readStatement(&is);
    } else {
        throw ParseError("WhileStatement::WhileStatement()",
                         "Invalid \"WHILE\" statement.");
    }
}


WhileStatement::~WhileStatement() {
    delete while_block;
}

#include "Utilities/OpalException.h"
void WhileStatement::execute(const Parser &parser) {
    curr = tokens.begin();
    keep = ++curr;
    Attribute condition = Attributes::makeBool("WHILE()", "while condition");

    try {
        condition.parse(*this, false);
        OpalData::getInstance()->update();

        while(Attributes::getBool(condition)) {
            while_block->execute(parser);
            OpalData::getInstance()->update();
        }
    } catch(...) {
        std::ostringstream oss;
        this->print(oss);
        ERRORMSG("Invalid WHILE condition '" + oss.str() + "'");
        throw;
    }
}