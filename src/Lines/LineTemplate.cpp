// ------------------------------------------------------------------------
// $RCSfile: LineTemplate.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LineTemplate
//   The class for storage of OPAL beam lines with arguments.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:41:39 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Lines/LineTemplate.h"
#include "AbstractObjects/OpalData.h"
#include "Lines/Line.h"
#include "Parser/SimpleStatement.h"
#include "Utilities/ParseError.h"
#include <vector>
#include <cassert>

// Class LineTemplate
// ------------------------------------------------------------------------

LineTemplate::LineTemplate():
    Macro(0, "LINE",
          "This object defines a beamline list with arguments.\n"
          "\t<name>(<args>);"),
    body("LINE")
{}


LineTemplate::LineTemplate(const std::string &name, Object *parent):
    Macro(name, parent), body(name)
{}


LineTemplate::~LineTemplate()
{}


LineTemplate *LineTemplate::clone(const std::string &/*name*/) {
    throw ParseError("LineTemplate::clone()",
                     "You cannot use this object without attributes.");
}


Object *LineTemplate::makeInstance
(const std::string &name, Statement &statement, const Parser *) {
    Line *instance = 0;

    try {
        // Parse actuals and check their number.
        parseActuals(statement);
        if(formals.size() != actuals.size()) {
            throw ParseError("LineTemplate::makeInstance()",
                             "Inconsistent number of macro arguments.");
        }

        // Expand the LINE macro in token form.
        body.start();
        Token token = body.readToken();
        SimpleStatement expansion(getOpalName(), 1);
        while(! token.isEOF()) {
            bool found = false;
            if(token.isWord()) {
                std::string word = token.getWord();

                for(std::vector<std::string>::size_type i = 0;
                    i < formals.size(); i++) {
                    if(word == formals[i]) {
                        std::vector<Token> act = actuals[i];
                        for(Token t : act) {
                            expansion.append(t);
                        }
                        found = true;
                        break;
                    }
                }
            }
            if(! found) expansion.append(token);
            token = body.readToken();
        }

        // Parse the modified line.
        Line *model = dynamic_cast<Line *>(OpalData::getInstance()->find("LINE"));
        instance = model->clone(name);
        instance->copyAttributes(*this);
        expansion.start();
        instance->parse(expansion);
    } catch(...) {
        delete instance;
        throw;
    }

    return instance;
}


Object *LineTemplate::makeTemplate(const std::string &, TokenStream &, Statement &) {
    // Should not be called.
    return 0;
}


void LineTemplate::parseTemplate(TokenStream &, Statement &statement) {
    parseFormals(statement);
    assert(statement.keyword("LINE"));

    // Store the template list.
    Token token = statement.getCurrent();
    if(token.isDel('=')) {
        body.append(token);
        int level = 0;
        while(! statement.atEnd()) {
            token = statement.getCurrent();
            if(token.isDel('(')) {
                level++;
            } else if(token.isDel(')')) {
                body.append(token);
                if(--level == 0) break;
            }
            body.append(token);
        }
    } else {
        throw ParseError("Line::makeTemplate()",
                         "Equals sign '=' expected.");
    }

}
