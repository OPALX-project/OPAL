#ifndef OPAL_RegularExpression_HH
#define OPAL_RegularExpression_HH 1

// ------------------------------------------------------------------------
// $RCSfile: RegularExpression.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RegularExpression.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <string>


// Class RegularExpression
// ------------------------------------------------------------------------
/// A regular expression.
//  This class encapsulates the UNIX regular expressions as defined in
//  regexp(5).  It provides a simple interface to the POSIX-compliant
//  regcomp/regexec package.

class RegularExpression {

public:

    /// Constructor.
    //  Construct regular expression from [b]pattern[/b].
    //  If [b]ignore[/b] is true, the expression ignores upper/lower case.
    RegularExpression(const std::string &pattern, bool ignore = false);

    RegularExpression(const RegularExpression &);
    ~RegularExpression();

    /// Match a string against the pattern.
    bool match(const std::string &s) const;

    /// Check the regular expression for sanity.
    bool OK() const;

private:

    // Not implemented.
    void operator = (const RegularExpression &);

    // Initialise.
    void init();

    // A copy of the pattern string.
    const std::string patt;

    // The flag to ignore case.
    bool caseIgnore;

    // The compiled regular expression.
    class Expression;
    Expression *expr;

    // The state flag, indicating any error condition.
    int state;
};

#endif // OPAL_RegularExpressions_HH
