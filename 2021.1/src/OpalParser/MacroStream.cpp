// ------------------------------------------------------------------------
// $RCSfile: MacroStream.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MacroStream
//   Implements an input buffer for reading tokens from a token list.
//   This string is intended for storing OPAL-9 macros.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/MacroStream.h"
#include "Parser/Token.h"


// Class MacroStream
// ------------------------------------------------------------------------

MacroStream::MacroStream(const std::string &macName):
    TokenStream(macName), body(), curr(body.begin())
{}


MacroStream::~MacroStream()
{}


void MacroStream::append(Token &token) {
    body.push_back(token);
}


Token MacroStream::readToken() {
    if(put_back_flag) {
        put_back_flag = false;
        return put_back;
    } else if(curr == body.end()) {
        return Token(stream_name, 1, Token::IS_EOF, "End of macro");
    } else {
        return *curr++;
    }
}


void MacroStream::start() {
    curr = body.begin();
    put_back_flag = false;
}
