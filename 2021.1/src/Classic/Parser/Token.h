#ifndef CLASSIC_Token_HH
#define CLASSIC_Token_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Token.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Token
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <iosfwd>
#include <string>


// Class Token
// ------------------------------------------------------------------------
/// Representation of a single input token.
//  All tokens contain the name of the input stream and the line number
//  where they came from.

class Token {

public:

    /// Possible token types.
    enum Type {
        IS_DELIMITER,
        IS_EOF,
        IS_ERROR,
        IS_INTEGER,
        IS_REAL,
        IS_WORD,
        IS_STRING
    };

    /// Constructor.
    //  Construct empty token.
    Token();

    /// Constructor.
    //  Construct character token with type [b]type[/b] and value [b]c[/b].
    Token(const std::string &file, int line, Type type, char c);

    /// Constructor.
    //  Construct std::string token with type [b]type[/b] and value [b]s[/b].
    Token(const std::string &file, int line, Type type, const char *s);

    /// Constructor.
    //  Construct string token with type [b]type[/b] and value [b]lex[/b].
    Token(const std::string &file, int line, Type type, const std::string &lex);

    /// Constructor.
    //  Construct real numeric token with lexeme [b]lex[/b] and value
    //  [b]value[/b].
    Token(const std::string &file, int line, const std::string &lex,
          double value);

    /// Constructor.
    //  Construct integer token with lexeme [b]lex[/b] and value [b]value[/b].
    Token(const std::string &file, int line, const std::string &lex, int value);

    Token(const Token &);
    ~Token();
    const Token &operator=(const Token &);

    /// Test for delimiter.
    //  Return true, if token is single character [b]del[/b].
    bool isDel(char del) const;

    /// Test for delimiter.
    //  Return true, if token is character string [b]del[/b].
    bool isDel(const char *del) const;

    /// Test for any delimiter.
    bool isDel()     const;

    /// Test for end of file.
    bool isEOF() const;

    /// Test for error.
    bool isError() const;

    /// Test for integer.
    bool isInteger() const;

    /// Test for real number.
    bool isReal() const;

    /// Test for word.
    bool isWord() const;

    /// Test for string.
    bool isString() const;

    /// Test for keyword.
    bool isKey(const char *key) const;

    /// Return boolean value.
    //  Throw ParseError, if token is not boolean.
    bool getBool() const;

    /// Return integer value.
    //  Throw ParseError, if token is not numeric.
    int getInteger() const;

    /// Return real value.
    //  Throw ParseError, if token is not numeric.
    double getReal() const;

    /// Return string value.
    //  Throw ParseError, if token is not string.
    std::string getString() const;

    /// Return word value.
    //  Throw ParseError, if token is not word.
    std::string getWord() const;

    /// Return the lexeme.
    const std::string &getLex() const;

    /// Return the token type.
    Type getType() const;

    /// Return the token's file name.
    const std::string &getFile() const;

    /// Return the token's line number.
    int getLine() const;

private:

    // Invalid token type.
    void invalid(const char *) const;

    // Input line and file.
    std::string file;
    int line;

    // Token type.
    Type type;

    // Lexeme for token.
    std::string lexeme;

    // Value for token.
    double   d_value;
    int      i_value;
    char     c_value;
};


// Output operator.
std::ostream &operator<<(std::ostream &, const Token &);

#endif // CLASSIC_Token_HH
