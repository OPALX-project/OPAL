#ifndef CLASSIC_AbsFileStream_HH
#define CLASSIC_AbsFileStream_HH

// ------------------------------------------------------------------------
// $RCSfile: AbsFileStream.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbsFileStream
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


// Class AbsFileStream
// ------------------------------------------------------------------------
/// A stream of input tokens.
//  The source of tokens is an abstract file stream.

class AbsFileStream: public TokenStream {

public:

    /// Constructor.
    //  Store stream name.
    explicit AbsFileStream(const std::string &name);

    virtual ~AbsFileStream();

    /// Read next input line.
    virtual bool fillLine() = 0;

    /// Read single token from file.
    //  Returns false when end of file is hit.
    virtual Token readToken();

protected:

    // Read double from current statement.
    Token readNumber();

    // Read string value from current statement.
    Token readString();

    // Read word from current statement.
    Token readWord();

    // Skip C-like comments.
    bool skipComments();

    // Current input line and position in line.
    std::string line;
    int curr_char;

private:

    // Not implemented.
    AbsFileStream();
    AbsFileStream(const AbsFileStream &);
    void operator=(const AbsFileStream &);
};

#endif // CLASSIC_AbsFileStream_HH