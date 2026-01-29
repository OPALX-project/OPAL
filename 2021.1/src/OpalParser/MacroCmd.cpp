//
// Class MacroCmd
//
//   This class parses the MACRO command.
//   Encapsulate the buffer for the ``archetypes'' of all macros.
//   The macro is stored as a MacroStream.  For execution, first the
//   parameters are replaced, then the resulting stream is sent to the parser.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "OpalParser/MacroCmd.h"

#include "Utility/PAssert.h"

#include "AbstractObjects/OpalData.h"
#include "OpalParser/OpalParser.h"
#include "Parser/Statement.h"
#include "Utilities/ParseError.h"
#include <vector>

MacroCmd::MacroCmd():
    Macro(0u, "MACRO",
          "A \"MACRO\" command defines a subroutine:\n"
          "\t<name>(<arguments>):MACRO{<body>}"),
    body(0)
{}


MacroCmd::MacroCmd(const std::string &name, MacroCmd *parent):
    Macro(name, parent), body() {
    body = new MacroStream(name);
}


MacroCmd::~MacroCmd()
{}


void MacroCmd::execute() {
    body->start();
    itsParser->run(&*body);
}


Object *MacroCmd::makeInstance
(const std::string &name, Statement &statement, const Parser *parser) {
    parseActuals(statement);

    // Check for consistency in argument number.
    if(formals.size() != actuals.size()) {
        throw ParseError("MacroCmd::makeInstance()",
                         "Inconsistent number of macro arguments.");
    }

    // Substitute the actual arguments.
    MacroCmd *macro = new MacroCmd(name, this);
    macro->itsParser = parser;
    body->start();
    Token token = body->readToken();

    while(! token.isEOF()) {
        bool found = false;

        if(token.isWord()) {
            std::string word = token.getWord();

            for(std::vector<std::string>::size_type i = 0;
                i < formals.size(); i++) {
                if(word == formals[i]) {
                    std::vector<Token> act = actuals[i];
                    for(Token t : act) {
                        macro->body->append(t);
                    }
                    found = true;
                    break;
                }
            }
        }

        if(! found) macro->body->append(token);
        token = body->readToken();
    }

    return macro;
}


Object *MacroCmd::makeTemplate
(const std::string &name, TokenStream &, Statement &statement) {
    MacroCmd *macro = new MacroCmd(name, this);
    macro->parseFormals(statement);

    // Parse macro body->
    bool isMacro = statement.keyword("MACRO");
    PAssert(isMacro);
    Token token;

    if(statement.delimiter('{')) {
        int level = 1;
        while(true) {
            if(statement.atEnd()) {
                throw ParseError("MacroCmd::makeTemplate()",
                                 "MACRO body is not closed.");
            } else {
                token = statement.getCurrent();

                if(token.isDel('{')) {
                    ++level;
                } else if(token.isDel('}')) {
                    if(--level == 0) break;
                }

                macro->body->append(token);
            }
        }
    } else {
        throw ParseError("MacroCmd::makeTemplate()",
                         "Missing MACRO body, should be \"{...}\".");
    }

    return macro;
}