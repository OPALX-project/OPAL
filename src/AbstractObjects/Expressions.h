//
// Namespace Expressions
//   This namespace contains the representations for all OPAL expressions
//   (scalar and array). It also defines the entire user interface for
//   the expression parser. All classes derived from these classes are
//   declared in this same namespace in the module "Expressions".
//
//   The expression parser is a recursive-descent parser.  When a parse
//   routine recognizes an expression, it returns its internal represenation
//   and skips the corresponding input tokens.  In case of error, a
//   ParseError exception is thrown.
//
// Template class: Expression::Scalar<T>
//   This abstract base class defines the interface for OPAL scalar
//   expressions whose values have type T.
//
// Template class: Expression::Array<T>
//   This abstract base class defines the interface for OPAL array
//   expressions whose components have type T.
//
// Template class: Expression::PtrToScalar<T>
//   This class implements a pointer which owns a scalar expression
//   of type Expression::Scalar<T>.
//
// Template class: Expression::PtrToArray<T>
//   This class implements a pointer which owns an array expression
//   of type Expression::Array<T>.
//
// Template class: Expression::ArrayOfPtrs<T>
//   This class implements an array of pointers, each of which owns
//   an expression of type Expression::Scalar<T>.
//
// Declared functions:
//   A collection of functions implementing a recursive descent parser.
//
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#ifndef OPAL_Expressions_HH
#define OPAL_Expressions_HH 1

#include "MemoryManagement/OwnPtr.h"
#include "Parser/Token.h"
#include <list>
#include <iosfwd>
#include <string>
#include <vector>

class PlaceRep;
class RangeRep;
class SFunction;
class Statement;
class Table;
class TableRowRep;

// ========================================================================
/// Representation objects and parsers for attribute expressions.

namespace Expressions {


    // Template class Expression::Scalar<T>.
    // ----------------------------------------------------------------------
    /// A scalar expression.

    template <class T> class Scalar {

    public:

        Scalar();
        Scalar(const Scalar &);
        virtual ~Scalar();

        /// Copy scalar expression.
        //  Deep copy of the expression.
        virtual Scalar<T> *clone() const = 0;

        /// Evaluate.
        //  Recursively evaluate the expression and return the result.
        //  Any non-deferred values are cached.
        virtual T evaluate() const = 0;

        /// Test for constant.
        //  Default action assumes non-constant.
        virtual bool isConstant() const;

        /// Print expression.
        //  Print the expression on the specified stream.
        //  The [b]precedence[/b] parameter controls printing of parentheses.
        virtual void print(std::ostream &, int precedence = 99) const = 0;

    private:

        // Not implemented.
        void operator=(const Scalar &);
    };


    // Template class Expression::PtrToScalar<T>.
    // ----------------------------------------------------------------------
    /// A pointer to a scalar expression.

    template <class T> class PtrToScalar:
        public OwnPtr<Scalar<T> > {

    public:

        /// Constructor from an object just created.
        PtrToScalar(Scalar<T> *rhs);

        PtrToScalar();
        PtrToScalar(const PtrToScalar &rhs);
        ~PtrToScalar();
        PtrToScalar& operator=(const PtrToScalar&) = default;
    };


    // Template class Expression::ArrayOfPtrs<T>.
    // ------------------------------------------------------------------------
    /// An array of pointers to scalar expressions.

    template <class T> class ArrayOfPtrs: public std::vector<PtrToScalar<T> > {

    public:

        ArrayOfPtrs();
        ArrayOfPtrs(int);
        ArrayOfPtrs(const std::vector<PtrToScalar<T> > &rhs);
        ~ArrayOfPtrs();
    };


    // Template class Expression::Array<T>.
    // ----------------------------------------------------------------------
    /// An array expression.

    template <class T> class OArray {

    public:

        OArray();
        OArray(const OArray &);
        virtual ~OArray();

        /// Copy expression.
        //  Deep copy.
        virtual OArray<T> *clone() const = 0;

        /// Evaluate.
        //  Recursively evaluate the expression and return the result.
        //  Any non-deferred expression is cached.
        virtual std::vector<T> evaluate() const = 0;

        /// Test for constant.
        //  Default action assumes non-constant.
        virtual bool isConstant() const;

        /// Print expression.
        //  Print the expression on the specified stream.
        //  The [b]precedence[/b] parameter controls printing of parentheses.
        virtual void print(std::ostream &, int precedence = 99) const = 0;

    private:

        // Not implemented.
        void operator=(const OArray &);
    };


    // Template class Expression::PtrToArray<T>.
    // ----------------------------------------------------------------------
    /// A pointer to an array expression.

    template <class T> class PtrToArray:
        public OwnPtr<OArray<T> > {

    public:

        /// Constructor from object just created.
        PtrToArray(OArray<T> *rhs);

        PtrToArray();
        PtrToArray(const PtrToArray &rhs);
        ~PtrToArray();
        PtrToArray& operator=(const PtrToArray<T>&) = default;
    };


    // ----------------------------------------------------------------------
    // SCALAR EXPRESSION PARSERS.

    /// Parse boolean expression.
    extern PtrToScalar<bool> parseBool(Statement &);

    /// Parse real expression.
    extern PtrToScalar<double> parseReal(Statement &);

    /// Parse real constant.
    extern double parseRealConst(Statement &);

    /// Parse string value.
    //  When no string is seen, a ParseError is thrown with the message
    //  given as the second argument.
    extern std::string parseString(Statement &, const char msg[]);


    // ARRAY EXPRESSION PARSERS.

    /// Parse boolean array expression.
    extern PtrToArray<bool> parseBoolArray(Statement &);

    /// Parse real array expression.
    extern PtrToArray<double> parseRealArray(Statement &);

    /// Parse real array constant.
    extern PtrToArray<double> parseRealConstArray(Statement &);

    /// Parse string array.
    extern std::vector<std::string> parseStringArray(Statement &);


    // SPECIAL VALUE PARSERS.

    /// Test for one-character delimiter.
    extern void parseDelimiter(Statement &stat, char delim);

    /// Test for two-character delimiter.
    extern void parseDelimiter(Statement &stat, const char delim[2]);

    /// Parse place specification.
    extern PlaceRep parsePlace(Statement &);

    /// Parse range specification.
    extern RangeRep parseRange(Statement &);

    // Forward declaration.
    template <class T> class SRefAttr;

    /// Parse variable reference.
    extern SRefAttr<double> *parseReference(Statement &);

    /// Parse a token list (for macro argument and the like).
    extern TableRowRep parseTableRow(Statement &);

    /// Parse a token list (for macro argument and the like).
    extern std::list<Token> parseTokenList(Statement &);

    /// Parse a token list array (for LIST commands).
    extern std::vector<std::list<Token> > parseTokenListArray(Statement &);


    // SPECIAL PARSER.

    /// Parse table expression (depends on a table's rows).
    extern PtrToScalar<double> parseTableExpression(Statement &, const Table *);


    // Implementation of class Expression::Scalar<T>.
    // ----------------------------------------------------------------------

    template <class T> inline
    Scalar<T>::Scalar()
    {}


    template <class T> inline
    Scalar<T>::Scalar(const Scalar<T> &)
    {}


    template <class T> inline
    Scalar<T>::~Scalar()
    {}


    template <class T> inline
    bool Scalar<T>::isConstant() const {
        return false;
    }


    // Implementation of class Expression::OArray<T>.
    // ----------------------------------------------------------------------

    template <class T> inline
    OArray<T>::OArray()
    {}


    template <class T> inline
    OArray<T>::OArray(const OArray<T> &)
    {}


    template <class T> inline
    OArray<T>::~OArray()
    {}


    template <class T> inline
    bool OArray<T>::isConstant() const {
        return false;
    }


    // Implementation of class Expression::PtrToScalar<T>.
    // ----------------------------------------------------------------------


    template <class T> inline
    PtrToScalar<T>::PtrToScalar():
        OwnPtr<Scalar<T> >()
    {}


    template <class T> inline
    PtrToScalar<T>::PtrToScalar(const PtrToScalar &rhs):
        OwnPtr<Scalar<T> >(rhs)
    {}


    template <class T> inline
    PtrToScalar<T>::PtrToScalar(Scalar<T> *rhs):
        OwnPtr<Scalar<T> >(rhs)
    {}


    template <class T> inline
    PtrToScalar<T>::~PtrToScalar()
    {}


    // Implementation of class Expression::ArrayOfPtrs<T>.
    // ----------------------------------------------------------------------

    template <class T> inline
    ArrayOfPtrs<T>::ArrayOfPtrs():
        std::vector<PtrToScalar<T> >()
    {}


    template <class T> inline
    ArrayOfPtrs<T>::ArrayOfPtrs(int n):
        std::vector<PtrToScalar<T> >(n, PtrToScalar<T>())
    {}


    template <class T> inline
    ArrayOfPtrs<T>::ArrayOfPtrs(const std::vector<PtrToScalar<T> > &rhs):
        std::vector<PtrToScalar<T> >(rhs)
    {}


    template <class T> inline
    ArrayOfPtrs<T>::~ArrayOfPtrs()
    {}


    // Implementation of class Expression::PtrToArray<T>.
    // ----------------------------------------------------------------------


    template <class T> inline
    PtrToArray<T>::PtrToArray():
        OwnPtr<OArray<T> >()
    {}


    template <class T> inline
    PtrToArray<T>::PtrToArray(const PtrToArray &rhs):
        OwnPtr<OArray<T> >(rhs)
    {}


    template <class T> inline
    PtrToArray<T>::PtrToArray(OArray<T> *rhs):
        OwnPtr<OArray<T> >(rhs)
    {}


    template <class T> inline
    PtrToArray<T>::~PtrToArray()
    {}


    // End of namespace Expressions.
    // ======================================================================
};

#endif // OPAL_Expressions_HH

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
