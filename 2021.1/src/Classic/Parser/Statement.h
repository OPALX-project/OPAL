#ifndef MAD_Statement_HH
#define MAD_Statement_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Statement.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Statement
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2002/10/05 00:48:29 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Parser/Token.h"
#include <iosfwd>
#include <list>
#include <string>

#include "Utility/Inform.h"

class Parser;


// class Statement
// ------------------------------------------------------------------------
/// Interface for statements.
//  The statement is stored as a list of Token's.

class Statement {

public:

    /// The type of the enclosed token list.
    typedef std::list<Token> TokenList;

    /// Constructor.
    //  Store the stream name and the line where the statement begins.
    Statement(const std::string &name, int line);

    /// Constructor.
    //  Stores a name (e.g. for a macro) and the token list.
    Statement(const std::string &name, TokenList &);

    virtual ~Statement();

    /// Append a token.
    void append(const Token &);

    /// Test for end of command.
    //  This method is called by the parser in order to find out wether
    //  the end of the statement has been reached.
    bool atEnd() const;

    /// Return boolean value.
    //  If the next item is a boolean literal:
    //  [ol][li]Set [b]value[/b] to its value.
    //  [li]Skip the value in the token list.
    //  [li]Return true.
    //  [/ol]
    //  Otherwise return false.
    bool boolean(bool &value);

    /// Test for delimiter.
    //  If the next item is the given character,
    //  skip it and return true, otherwise return false.
    bool delimiter(char c);

    /// Test for delimiter choice.
    //  If the next item is one of the characters in the given string,
    //  skip it and return true, otherwise return false.
    bool delimiter(const char *s);

    /// Execute.
    //  This method must be specially defined in conditional or loop
    //  statements. Normally it just calls the parser to execute the
    //  statement.
    virtual void execute(const Parser &) = 0;

    /// Return current token and skip it.
    Token &getCurrent();

    /// Return signed integer.
    //  If the next item is an integer:
    //  [ol][li]Set [b]value[/b] to its value.
    //  [li]Skip the value in the token list.
    //  [li]Return true.
    //  [/ol]
    //  Otherwise return false.
    bool integer(int &value);

    /// Return unsigned integer.
    //  If the next item is an integer:
    //  [ol][li]Set [b]value[/b] to its value.
    //  [li]Skip the value in the token list.
    //  [li]Return true.
    //  [/ol]
    //  Otherwise return false.
    bool integer(unsigned &value);

    /// Test for keyword.
    //  If the next item is the keyword [b]s[/b], skip it and return true,
    //  otherwise return false.
    bool keyword(const char *s);

    /// Return real value.
    //  If the next item is a real number:
    //  [ol][li]Set [b]value[/b] to its value.
    //  [li]Skip the value in the token list.
    //  [li]Return true.
    //  [/ol]
    //  Otherwise return false.
    bool real(double &value);

    /// Return string value.
    //  If the next item is a string literal:
    //  [ol][li]Set [b]value[/b] to its value.
    //  [li]Skip the value in the token list.
    //  [li]Return true.
    //  [/ol]
    //  Otherwise return false.
    bool str(std::string &value);

    /// Return word value.
    //  If the next item is a word:
    //  [ol][li]Set [b]value[/b] to its value.
    //  [li]Skip the value in the token list.
    //  [li]Return true.
    //  [/ol]
    //  Otherwise return false.
    bool word(std::string &value);

    /// Mark position in command.
    //  Parsing can later be resumed by calling restore().
    void mark();

    /// Return to marked position.
    //  Resume parsing at the position wher mark() was last called.
    void restore();

    /// Return to start.
    //  Resume parsing at the beginning of the statement.
    void start();

    /// Skip.
    //  Skip tokens up to next comma or end of statement, whichever comes first.
    void skip();

    /// Return current character number in line
    unsigned int position() const;

    /// Print statement.
    //  Print the statement on [b]os[/b].
    virtual void print(std::ostream &os) const;

    /// Print position.
    //  Print a message, containing the stream name and its line in the input
    //  stream.  If [b]withToken[/b] is true, print also the last token parsed.
    virtual void printWhere(Inform &msg, bool withToken) const;

    std::string str() const;
protected:

    // Line number where statement begins.
    int stat_line;

    // Input stream name.
    std::string buffer_name;

    // Token list.
    TokenList tokens;
    TokenList::iterator curr;
    TokenList::iterator keep;
};


// Output operator.
inline std::ostream &operator<<(std::ostream &os, const Statement &statement) {
    statement.print(os);
    return os;
}

inline Inform &operator<<(Inform &os, const Statement &statement) {
    std::ostringstream msg;
    statement.print(msg);
    os << msg.str();

    return os;
}

#endif // MAD_Statement_HH