// ------------------------------------------------------------------------
// $RCSfile: SequenceTemplate.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SequenceTemplate
//   The class for storage of OPAL beam sequences with arguments.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:41:41 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Lines/SequenceTemplate.h"
#include "AbstractObjects/OpalData.h"
#include "Lines/Sequence.h"
#include "Lines/SequenceParser.h"
#include "OpalParser/MacroStream.h"
#include "OpalParser/OpalParser.h"
#include "Parser/SimpleStatement.h"
#include "Utilities/ParseError.h"
#include <vector>
#include <cassert>

// Class SequenceTemplate
// ------------------------------------------------------------------------

SequenceTemplate::SequenceTemplate():
    Macro(0, "SEQUENCE",
          "This object defines a beamsequence list with arguments.\n"
          "\t<name>(<args>);"),
    body("SEQUENCE")
{}


SequenceTemplate::SequenceTemplate(const std::string &name, Object *parent):
    Macro(name, parent), body(name)
{}


SequenceTemplate::~SequenceTemplate()
{}


SequenceTemplate *SequenceTemplate::clone(const std::string &/*name*/) {
    throw ParseError("SequenceTemplate::clone()",
                     "You cannot use this object without attributes.");
}


Object *SequenceTemplate::makeInstance
(const std::string &name, Statement &statement, const Parser *) {
    MacroStream *expansion = 0;
    Sequence *instance = 0;

    try {
        // Parse actuals and check their number.
        parseActuals(statement);
        if(formals.size() != actuals.size()) {
            throw ParseError("MacroCmd::makeInstance()",
                             "Inconsistent number of macro arguments.");
        }

        // Expand the SEQUENCE macro in token form.
        body.start();
        Token token = body.readToken();
        expansion = new MacroStream(getOpalName());
        while(! token.isEOF()) {
            bool found = false;
            if(token.isWord()) {
                std::string word = token.getWord();
                for(std::vector<std::string>::size_type i = 0;
                    i < formals.size(); i++) {
                    if(word == formals[i]) {
                        std::vector<Token> act = actuals[i];
                        for(Token t : act) {
                            expansion->append(t);
                        }
                        found = true;
                        break;
                    }
                }
            }
            if(! found) expansion->append(token);
            token = body.readToken();
        }

        // Make the instance and parse it.
        Sequence *model = dynamic_cast<Sequence *>(OpalData::getInstance()->find("SEQUENCE"));
        instance = model->clone(name);
        instance->copyAttributes(*this);
        expansion->start();
        SequenceParser parser(instance);
        parser.run(&*expansion);
    } catch(...) {
        delete expansion;
        delete instance;
        throw;
    }

    return instance;
}


Object *SequenceTemplate::makeTemplate
(const std::string &, TokenStream &, Statement &) {
    // Should not be called.
    return 0;
}


void SequenceTemplate::parseTemplate(TokenStream &is, Statement &statement) {
    // Save the formals.
    parseFormals(statement);
    assert(statement.keyword("SEQUENCE"));

    // Parse the sequence header.
    Object::parse(statement);

    // Save the sequence body.
    Token token = is.readToken();

    // Read through ENDSEQUENCE.
    while(! token.isEOF()) {
        body.append(token);
        if(token.isKey("ENDSEQUENCE")) {
            // Read remainder up to ';'.
            token = is.readToken();
            while(! token.isEOF()) {
                body.append(token);
                if(token.isDel(';')) break;
                token = is.readToken();
            }
            break;
        }
        token = is.readToken();
    }
}
