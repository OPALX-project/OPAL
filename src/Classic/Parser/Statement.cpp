// ------------------------------------------------------------------------
// $RCSfile: Statement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Statement
//   An abstract base class for all input language statements.
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
#include "Parser/Token.h"
#include <iostream>

// Class Statement
// ------------------------------------------------------------------------

Statement::Statement(const std::string &name, int line):
    stat_line(line), buffer_name(name), tokens()
{}


Statement::Statement(const std::string &name, TokenList &list):
    stat_line(1), buffer_name(name), tokens() {
    tokens.swap(list);
    curr = tokens.begin();
}


Statement::~Statement() {
    tokens.erase(tokens.begin(), tokens.end());
}


void Statement::append(const Token &token) {
    tokens.push_back(token);
}


bool Statement::atEnd() const {
    return TokenList::const_iterator(curr) == tokens.end();
}


bool Statement::boolean(bool &value) {
    if(curr != tokens.end()  &&  curr->isWord()) {
        std::string word = curr->getWord();

        if(word == "TRUE") {
            value = true;
            ++curr;
            return true;
        } else if(word == "FALSE") {
            value = false;
            ++curr;
            return true;
        }
    }

    return false;
}


Token &Statement::getCurrent() {
    return *(curr++);
}


bool Statement::integer(int &value) {
    if(curr != tokens.end()  &&  curr->isInteger()) {
        value = curr->getInteger();
        ++curr;
        return true;
    } else {
        return false;
    }
}


bool Statement::integer(unsigned &value) {
    if(curr != tokens.end()  &&  curr->isInteger()) {
        value = curr->getInteger();
        ++curr;
        return true;
    } else {
        return false;
    }
}


bool Statement::delimiter(char c) {
    if(curr != tokens.end()  && (*curr).isDel(c)) {
        ++curr;
        return true;
    } else {
        return false;
    }
}


bool Statement::delimiter(const char *s) {
    if(curr != tokens.end()  && (*curr).isDel(s)) {
        ++curr;
        return true;
    } else {
        return false;
    }
}


bool Statement::keyword(const char *key) {
    if(curr != tokens.end()  && (*curr).isKey(key)) {
        ++curr;
        return true;
    } else {
        return false;
    }
}


bool Statement::real(double &value) {
    if(curr != tokens.end()) {
        if(curr->isReal()) {
            value = curr->getReal();
            ++curr;
            return true;
        } else if(curr->isInteger()) {
            value = double(curr->getInteger());
            ++curr;
            return true;
        }
    }

    return false;
}


bool Statement::str(std::string &value) {
    if(curr != tokens.end()  &&  curr->isString()) {
        value = curr->getLex();
        ++curr;
        return true;
    } else {
        return false;
    }
}


bool Statement::word(std::string &value) {
    if(curr != tokens.end()  &&  curr->isWord()) {
        value = curr->getLex();
        ++curr;
        return true;
    } else {
        return false;
    }
}


void Statement::mark() {
    keep = curr;
}


void Statement::restore() {
    curr = keep;
}


void Statement::start() {
    curr = tokens.begin();
}


void Statement::skip() {
    while(! atEnd()  &&  !(*curr).isDel(','))
        curr++;
}

unsigned int Statement::position() const {
    std::ostringstream os;
    bool white = false;

    for(TokenList::const_iterator c = tokens.begin(); c != curr; ++c) {
        if(white && !c->isDel()) os << ' ';
        white = !c->isDel();
        os << *c;
    }
    if(white && !std::next(curr)->isDel()) os << ' ';

    return os.str().length() - 1;
}

void Statement::print(std::ostream &msg) const {
    bool white = false;

    for(TokenList::const_iterator c = tokens.begin(); c != tokens.end(); ++c) {
        if(white && !c->isDel()) msg << ' ';
        white = !c->isDel();
        msg << *c;
    }

    msg << ';' << std::endl;
}


void Statement::printWhere(Inform &msg, bool withToken) const {

    msg << "*** in line " << stat_line << " of file \"" << buffer_name << "\"";

    if(withToken) {
        if(TokenList::const_iterator(curr) == tokens.end()) {
            msg << " at end of statement:" << endl;
        } else {
            msg << " before token \"" << *curr << "\":" << endl;
        }
    } else {
        msg << ":\n";
    }
}

std::string Statement::str() const {
    std::ostringstream str;
    print(str);

    return str.str();
}