// ------------------------------------------------------------------------
// $RCSfile: TokenStream.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TokenStream
//   Abstract base class for input token streams.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Parser/TokenStream.h"
#include "Utilities/LogicalError.h"
#include <cctype>


// Class TokenStream
// ------------------------------------------------------------------------

TokenStream::TokenStream(const std::string &name):
    stream_name(name), curr_line(0), put_back_flag(false), put_back()
{}


TokenStream::~TokenStream()
{}


void TokenStream::putBack(const Token &token) {
    if(put_back_flag) {
        throw LogicalError("TokenStream::pushBack()",
                           "Cannot push back two tokens.");
    }

    put_back = token;
    put_back_flag = true;
}


int TokenStream::getLine() const {
    return curr_line;
}


const std::string &TokenStream::getName() const {
    return stream_name;
}
