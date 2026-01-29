// ------------------------------------------------------------------------
// $RCSfile: Token.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Token:
//   Holds a single input token.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Parser/Token.h"
#include "Utilities/ParseError.h"
#include <cctype>
#include <iostream>
#include <sstream>


// Class Token
// ------------------------------------------------------------------------

Token::Token():
    file(), line(0), type(IS_ERROR), lexeme(),
    d_value(0.0), i_value(0), c_value(0)
{}


Token::Token(const Token &rhs):
    file(rhs.file), line(rhs.line), type(rhs.type), lexeme(rhs.lexeme),
    d_value(rhs.d_value), i_value(rhs.i_value), c_value(rhs.c_value)
{}


Token::Token(const std::string &fil, int lin, Type typ, char value):
    file(fil), line(lin), type(typ), lexeme(1u, value),
    d_value(0.0), i_value(0), c_value(value)
{}


Token::Token(const std::string &fil, int lin, Type typ, const char *value):
    file(fil), line(lin), type(typ), lexeme(value),
    d_value(0.0), i_value(0), c_value(0)
{}


Token::Token
(const std::string &fil, int lin, Type typ, const std::string &lex):
    file(fil), line(lin), type(typ), lexeme(lex),
    d_value(0.0), i_value(0), c_value(0)
{}


Token::Token
(const std::string &fil, int lin, const std::string &lex, double value):
    file(fil), line(lin), type(IS_REAL), lexeme(lex),
    d_value(value), i_value(0), c_value(0)
{}


Token::Token
(const std::string &fil, int lin, const std::string &lex, int value):
    file(fil), line(lin), type(IS_INTEGER), lexeme(lex),
    d_value(0.0), i_value(value), c_value(0)
{}


Token::~Token()
{}


const Token &Token::operator=(const Token &rhs) {
    file = rhs.file;
    line = rhs.line;
    type = rhs.type;
    lexeme = rhs.lexeme;
    d_value = rhs.d_value;
    i_value = rhs.i_value;
    c_value = rhs.c_value;
    return *this;
}


bool Token::isDel(char del) const {
    return (type == IS_DELIMITER  &&  lexeme == std::string(1u, del));
}


bool Token::isDel(const char *del) const {
    return (type == IS_DELIMITER  &&  lexeme == del);
}


bool Token::isDel() const {
    return (type == IS_DELIMITER);
}


bool Token::isEOF() const {
    return (type == IS_EOF);
}


bool Token::isError() const {
    return (type == IS_ERROR);
}


bool Token::isInteger() const {
    return (type == IS_INTEGER);
}


bool Token::isReal() const {
    return (type == IS_REAL);
}


bool Token::isWord() const {
    return (type == IS_WORD);
}


bool Token::isString() const {
    return (type == IS_STRING);
}


bool Token::isKey(const char *key) const {
    return (type == IS_WORD  &&  lexeme == key);
}


bool Token::getBool() const {
    if(type == IS_WORD) {
        if(lexeme == "TRUE") {
            return true;
        } else if(lexeme == "FALSE") {
            return false;
        }
    }

    invalid("boolean");
    return false;
}


int Token::getInteger() const {
    if(type == IS_INTEGER) {
        return i_value;
    } else if(type == IS_REAL) {
        return int(d_value + 0.5);
    }

    invalid("integer");
    return 0;
}


double Token::getReal() const {
    if(type == IS_INTEGER) {
        return double(i_value);
    } else if(type == IS_REAL) {
        return d_value;
    }

    invalid("real");
    return 0.0;
}


std::string Token::getString() const {
    if(type == IS_STRING) {
        return lexeme;
    }

    invalid("string");
    return "";
}


std::string Token::getWord() const {
    if(type == IS_WORD) {
        return lexeme;
    }

    invalid("identifier");
    return "";
}


const std::string &Token::getLex() const {
    return lexeme;
}


Token::Type Token::getType() const {
    return type;
}


const std::string &Token::getFile() const {
    return file;
}


int Token::getLine() const {
    return line;
}


std::ostream &operator<<(std::ostream &os, const Token &token) {
    switch(token.getType()) {

        case Token::IS_DELIMITER:
            os << token.getLex();
            break;

        case Token::IS_EOF:
            os << "End of file";
            break;

        case Token::IS_ERROR:
            os << "Error token: " << token.getLex();
            break;

        case Token::IS_INTEGER:
            os << token.getInteger();
            break;

        case Token::IS_REAL:
            os << token.getReal();
            break;

        case Token::IS_WORD:
            os << token.getLex();
            break;

        case Token::IS_STRING:
            os << '"' << token.getLex() << '"';
            break;
    }

    return os;
}


void Token::invalid(const char *type) const {
    std::ostringstream os;
    os << "Expected " << type << " token in line " << line << " of file \""
       << file << "\".";
    throw ParseError("Token::invalid()", os.str().c_str());
}
