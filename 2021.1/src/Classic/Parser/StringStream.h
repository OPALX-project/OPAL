#ifndef CLASSIC_StringStream_HH
#define CLASSIC_StringStream_HH

// ------------------------------------------------------------------------
// $RCSfile: StringStream.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StringStream
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


// Class StringStream
// ------------------------------------------------------------------------
/// A stream of input tokens.
//  The source of tokens is a C++ string.

class StringStream: public TokenStream {

public:

    /// Constructor.
    //  Initialize stream to read from the named file.
    StringStream(const std::string &s);

    virtual ~StringStream();

    /// Read single token from file.
    virtual Token readToken();

private:

    // Not implemented.
    StringStream();
    StringStream(const StringStream &);
    void operator=(const StringStream &);

    // Read double from current statement.
    Token readNumber();

    // Read string value from current statement.
    Token readString();

    // Read word from current statement.
    Token readWord();

    // Current input string;
    const std::string line_m;
    std::string::size_type currentChar_m;
};

#endif // CLASSIC_StringStream_HH