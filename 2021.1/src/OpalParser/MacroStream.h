#ifndef CLASSIC_MacroStream_HH
#define CLASSIC_MacroStream_HH

// ------------------------------------------------------------------------
// $RCSfile: MacroStream.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MacroStream
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Parser/TokenStream.h"
#include "Parser/Token.h"
#include <list>


// Class MacroStream
// ------------------------------------------------------------------------
/// An input buffer for macro commands.
//  This list of tokens stores a macro template or instantiation.

class MacroStream: public TokenStream {

public:

    /// Constructor.
    //  Assign the macro name as a buffer name.
    MacroStream(const std::string &);

    virtual ~MacroStream();

    /// Append a token to the stream.
    void append(Token &);

    /// Read a token from the stream.
    virtual Token readToken();

    /// Reset stream to start.
    void start();

private:

    // Not implemented.
    MacroStream();
    MacroStream(const MacroStream &);
    void operator=(const MacroStream &);

    // The list of tokens.
    typedef std::list<Token> TokenList;
    TokenList body;
    TokenList::iterator curr;
};

#endif // CLASSIC_MacroStream_HH
