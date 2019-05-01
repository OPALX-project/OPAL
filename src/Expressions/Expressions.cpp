// ------------------------------------------------------------------------
// $RCSfile: Expressions.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace Expressions:
//   This namespace contains the attribute parsers for OPAL:
//   a collection of functions.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algorithms/PartBunchBase.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/RangeRep.h"
#include "AbstractObjects/Table.h"
#include "AbstractObjects/TableRowRep.h"
#include "AbstractObjects/OpalData.h"
#include "Expressions/ABinary.h"
#include "Expressions/AColumn.h"
#include "Expressions/AList.h"
#include "Expressions/ARefExpr.h"
#include "Expressions/ARow.h"
#include "Expressions/ASUnary.h"
#include "Expressions/AUnary.h"
#include "Expressions/ATable.h"
#include "Expressions/Indexer.h"
#include "Expressions/SBinary.h"
#include "Expressions/SCell.h"
#include "Expressions/SConstant.h"
#include "Expressions/SFunction.h"
#include "Expressions/SHash.h"
#include "Expressions/SNull.h"
#include "Expressions/SRefAttr.h"
#include "Expressions/SRefExpr.h"
#include "Expressions/SUnary.h"
#include "Expressions/TFind.h"
#include "Expressions/TFunction0.h"
#include "Expressions/TFunction1.h"
#include "Expressions/TFunction2.h"
#include "Parser/Statement.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include "Utilities/CLRangeError.h"
#include "Utilities/Round.h"
#include "Utilities/Truncate.h"
#include "ValueDefinitions/BoolConstant.h"
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <list>
#include <sstream>
#include <vector>


// Namespace Expressions
// ------------------------------------------------------------------------

namespace Expressions {


    // Boolean operators.
    // ----------------------------------------------------------------------

    bool Or(bool a, bool b) { return a || b; }
    static const TFunction2<bool, bool> logOr = { "||", 1, Or };

    bool And(bool a, bool b) { return a && b; }
    static const TFunction2<bool, bool> logAnd = { "&&", 2, And };

    bool Le(double a, double b) { return a <= b; }
    static const TFunction2<bool, double> lessEqual = { "<=", 3, Le } ;

    bool Lt(double a, double b) { return a < b; }
    static const TFunction2<bool, double> less = { "<", 3, Lt } ;

    bool Ge(double a, double b) { return a >= b; }
    static const TFunction2<bool, double> greaterEqual = { ">=", 3, Ge } ;

    bool Gt(double a, double b) { return a > b; }
    static const TFunction2<bool, double> greater = { ">", 3, Gt } ;

    bool Eq(double a, double b) { return a == b; }
    static const TFunction2<bool, double> equal = { "==", 3, Eq } ;

    bool Ne(double a, double b) { return a != b; }
    static const TFunction2<bool, double> notEqual = { "!=", 3, Ne } ;

    // Real operators with one or two operands.
    // ----------------------------------------------------------------------

    double Neg(double a) {
        return - a;
    }

    static const TFunction1<double, double> negate = { "-", 5, Neg };

    double Sign(double a) {
        return (a > 0.) ? 1. : (a < 0.) ? -1. : 0.;
    }

    double Tgauss(double a) {
        double x, y;
        do { Options::rangen.gauss(x, y); }
        while(std::abs(x) > a);
        return x;
    }

    double Add(double a, double b) {
        return a + b;
    }

    static const TFunction2<double, double> plus =   { "+", 4, Add };

    double Sub(double a, double b) {
        return a - b;
    }

    static const TFunction2<double, double> minus =  { "-", 4, Sub };

    double Mpy(double a, double b) {
        return a * b;
    }

    static TFunction2<double, double> times = { "*", 5, Mpy };

    double Div(double a, double b) {
        if(b == 0.0) {
            errno = EDOM;
            return 1.0;
        } else {
            return a / b;
        }
    }

    static TFunction2<double, double> divide = { "/", 5, Div };

    // Real functions without arguments.
    // ----------------------------------------------------------------------


    double getEkin() {
        auto p = OpalData::getInstance()->getPartBunch();
        if (p)
            return p->get_meanKineticEnergy();
        else
            return -1.0;
    }


    double ranf() {
        return Options::rangen.uniform();
    }

    double gauss() {
        double x, y;
        Options::rangen.gauss(x, y);
        return x;
    }

    static const TFunction0<double> table0[] = {
        {"GETEKIN", -2, getEkin},
        {"RANF",    -2, ranf   },
        {"GAUSS",   -2, gauss  },
        {"SI",      -2, SFunction::arcIn  },
        {"SC",      -2, SFunction::arcCtr },
        {"SO",      -2, SFunction::arcOut },
        { 0,        -1, 0      }
    };

    // Real functions with two arguments.
    // ----------------------------------------------------------------------

    static const TFunction1<double, double> table1[] = {
        { "TRUNC",  -1, Truncate   },
        { "ROUND",  -1, Round      },
        { "FLOOR",  -1, std::floor },
        { "CEIL",   -1, std::ceil  },
        { "SIGN",   -1, Sign       },
        { "SQRT",   -1, std::sqrt  },
        { "LOG",    -1, std::log   },
        { "EXP",    -1, std::exp   },
        { "SIN",    -1, std::sin   },
        { "COS",    -1, std::cos   },
        { "ABS",    -1, std::abs   },
        { "TAN",    -1, std::tan   },
        { "ASIN",   -1, std::asin  },
        { "ACOS",   -1, std::acos  },
        { "ATAN",   -1, std::atan  },
        { "TGAUSS", -2, Tgauss     },
        { "ERF",    -1, std::erf   },
        { 0,        -1, 0          }
    };

    // Real functions with two arguments.
    // ----------------------------------------------------------------------

    double Max(double a, double b) {
        return std::max(a, b);
    }

    double Min(double a, double b) {
        return std::min(a, b);
    }

    double Mod(double a, double b) {
        if(b <= 0.0) errno = EDOM;
        return fmod(a, b);
    }

    static TFunction2<double, double> power = { "^", 6, pow };

    static const TFunction2<double, double> table2[] = {
        { "ATAN2",  -1, atan2  },
        { "MAX",    -1, Max    },
        { "MIN",    -1, Min    },
        { "MOD",    -1, Mod    },
        { "POW",    -2, pow    },
        { 0,        -1, 0      }
    };


    // Function of arrays giving a real result.
    // ----------------------------------------------------------------------

    double Mina(const std::vector<double> &array) {
        if(array.size()) {
            double result = array[0];
            for(std::vector<double>::size_type i = 1; i < array.size(); ++i) {
                result = Min(array[i], result);
            }
            return result;
        } else if(Options::warn) {
            std::cerr << "\n### Warning ### \"VMIN\" function of empty array.\n"
                      << std::endl;
        }
        return 0.0;
    }

    double Maxa(const std::vector<double> &array) {
        if(array.size()) {
            double result = array[0];
            for(std::vector<double>::size_type i = 1; i < array.size(); ++i) {
                result = std::max(array[i], result);
            }
            return result;
        } else if(Options::warn) {
            std::cerr << "\n### Warning ### \"VMAX\" function of empty array.\n"
                      << std::endl;
        }
        return 0.0;
    }

    double Rmsa(const std::vector<double> &array) {
        if(array.size()) {
            double result = 0.0;
            for(std::vector<double>::size_type i = 0; i < array.size(); ++i) {
                result += array[i] * array[i];
            }
            return sqrt(result / double(array.size()));
        } else if(Options::warn) {
            std::cerr << "\n### Warning ### \"VRMS\" function of empty array.\n"
                      << std::endl;
        }
        return 0.0;
    }

    double AbsMax(const std::vector<double> &array) {
        if(array.size()) {
            double result = std::abs(array[0]);
            for(std::vector<double>::size_type i = 1; i < array.size(); ++i) {
                result = std::max(std::abs(array[i]), result);
            }
            return result;
        } else if(Options::warn) {
            std::cerr << "\n### Warning ### \"VABSMAX\" function of empty array.\n"
                      << std::endl;
        }
        return 0.0;
    }

    typedef TFunction1<double, const std::vector<double>&> ArrayFun;
    static const ArrayFun tablea[] = {
        { "VMIN",    -1, Mina   },
        { "VMAX",    -1, Maxa   },
        { "VRMS",    -1, Rmsa   },
        { "VABSMAX", -1, AbsMax },
        { 0,         -1, 0      }
    };


    // Internal variables.
    const Table *currentTable = 0;
    OwnPtr<ATable> currentArray;


    // FORWARD DECLARATIONS FOR LOCAL FUNCTIONS.
    // ------------------------------------------------------------------------

    // Parse an "&&" expression.
    PtrToScalar<bool> parseAnd(Statement &stat);

    // Parse an array factor.
    PtrToArray<double> parseArrayFactor(Statement &stat);

    // Parse an array primary.
    PtrToArray<double> parseArrayPrimary(Statement &stat);

    // Parse an array term.
    PtrToArray<double> parseArrayTerm(Statement &stat);

    // Parse a token list in brackets (for macro argument and the like).
    // This method also checks the nesting of parentheses.
    void parseBracketList
    (Statement &stat, char close, std::list<Token> &result);

    // Parse a real factor.
    PtrToScalar<double> parseFactor(Statement &stat);

    // Parse a real primary expression.
    PtrToScalar<double> parsePrimary(Statement &stat);

    // Parse a relational expression.
    PtrToScalar<bool> parseRelation(Statement &stat);

    // Parse a real term.
    PtrToScalar<double> parseTerm(Statement &stat);

    // Parse a column generator "COLUMN(table,column,range)".
    PtrToArray<double> parseColumnGenerator(Statement &stat);

    // Parse a row generator "ROW(table,place,...)".
    PtrToArray<double> parseRowGenerator(Statement &stat);

    // Parse a table generator "TABLE(n1:n2:n3,expr)".
    PtrToArray<double> parseTableGenerator(Statement &stat);


    // PARSE SCALAR EXPRESSIONS.
    // ------------------------------------------------------------------------

    PtrToScalar<bool> parseBool(Statement &stat) {
        PtrToScalar<bool> second;
        PtrToScalar<bool> result = parseAnd(stat);

        while(stat.delimiter("||")) {
            second = parseAnd(stat);
            result = SBinary<bool, bool>::make(logOr, result, second);
        }

        return result;
    }


    PtrToScalar<double> parseReal(Statement &stat) {
        PtrToScalar<double> result;
        PtrToScalar<double> second;

        if(stat.delimiter('+')) {
            // Unary plus.
            result = parseTerm(stat);
        } else if(stat.delimiter('-')) {
            // Unary minus.
            result = parseTerm(stat);
            result = SUnary<double, double>::make(negate, result);
        } else {
            result = parseTerm(stat);
        }

        while(true) {
            if(stat.delimiter('+')) {
                second = parseTerm(stat);
                result = SBinary<double, double>::make(plus, result, second);
            } else if(stat.delimiter('-')) {
                second = parseTerm(stat);
                result = SBinary<double, double>::make(minus, result, second);
            } else {
                break;
            }
        }

        return result;
    }


    double parseRealConst(Statement &stat) {
        PtrToScalar<double> expr = parseReal(stat);
        return expr->evaluate();
    }


    std::string parseString(Statement &stat, const char msg[]) {
        std::string result = std::string("");
        std::string temp;

        while(stat.word(temp) || stat.str(temp)) {
            if(temp == "STRING") {
                parseDelimiter(stat, '(');
                double value = parseRealConst(stat);
                parseDelimiter(stat, ')');
                std::ostringstream os;
                os << value << std::ends;
                result += os.str();
            } else {
                result += temp;
            }

            if(! stat.delimiter('&')) break;
        }

        if(result.empty()) {
            throw ParseError("Expressions::parseString()", msg);
        } else {
            return result;
        }
    }


    // PARSE ARRAY EXPRESSIONS
    // ------------------------------------------------------------------------

    PtrToArray<bool> parseBoolArray(Statement &stat) {
        ArrayOfPtrs<bool> array;

        if(stat.delimiter('{')) {
            // List of boolean expressions within braces.
            do {
                array.push_back(parseBool(stat));
            } while(stat.delimiter(','));

            parseDelimiter(stat, '}');
        } else {
            // A single boolean expression.
            array.push_back(parseBool(stat));
        }

        return new AList<bool>(array);
    }


    PtrToArray<double> parseRealArray(Statement &stat) {
        PtrToArray<double> result;
        PtrToArray<double> second;

        if(stat.delimiter('+')) {
            // Unary plus.
            result = parseArrayTerm(stat);
        } else if(stat.delimiter('-')) {
            // Unary minus.
            result = parseArrayTerm(stat);
            result = new AUnary<double, double>(negate, result);
        } else {
            result = parseArrayTerm(stat);
        }

        while(true) {
            if(stat.delimiter('+')) {
                second = parseArrayTerm(stat);
                result = new ABinary<double, double>(plus, result, second);
            } else if(stat.delimiter('-')) {
                second = parseArrayTerm(stat);
                result = new ABinary<double, double>(minus, result, second);
            } else {
                break;
            }
        }

        return result;
    }


    PtrToArray<double> parseRealConstArray(Statement &stat) {
        PtrToArray<double> arg = parseRealArray(stat);
        std::vector<double> value = arg->evaluate();
        ArrayOfPtrs<double> array;

        for(ArrayOfPtrs<double>::size_type i = 0; i < value.size(); ++i) {
            array.push_back(new SConstant<double>(value[i]));
        }

        return new AList<double>(array);
    }


    std::vector<std::string> parseStringArray(Statement &stat) {
        std::vector<std::string> array;

        if(stat.delimiter('{')) {
            // List of string values within braces.
            do {
                array.push_back(parseString(stat, "String value expected."));
            } while(stat.delimiter(','));

            parseDelimiter(stat, '}');
        } else {
            // A single string value.
            array.push_back(parseString(stat, "String value expected."));
        }

        return array;
    }


    void parseDelimiter(Statement &stat, char delim) {
        if(! stat.delimiter(delim)) {
            throw ParseError("Expressions::parseDelimiter()",
                             std::string("Delimiter '") + delim + "' expected.");
        }
    }


    void parseDelimiter(Statement &stat, const char delim[2]) {
        if(! stat.delimiter(delim)) {
            throw ParseError("Expressions::parseDelimiter()",
                             std::string("Delimiter '") + delim + "' expected.");
        }
    }


    PlaceRep parsePlace(Statement &stat) {
        if(stat.keyword("SELECTED")) {
            return PlaceRep("SELECTED");
        } else {
            PlaceRep pos;

            do {
                if(stat.delimiter('#')) {
                    static char msg[] = "Expected 'S' or 'E' after '#'.";
                    if(stat.keyword("S")) {
                        pos.append("#S", 0);
                    } else if(stat.keyword("E")) {
                        pos.append("#E", 0);
                    } else {
                        throw ParseError("Expression::parsePlace()", msg);
                    }
                    break;
                } else {
                    std::string name = parseString(stat, "Expected <name> or '#'.");
                    int occurrence = 1;
                    if(stat.delimiter('[')) {
                        occurrence = int(Round(parseRealConst(stat)));
                        parseDelimiter(stat, ']');

                        if(occurrence <= 0) {
                            throw ParseError("Expressions::parsePlace()",
                                             "Occurrence counter must be positive.");
                        }
                    }
                    pos.append(name, occurrence);
                }
            } while(stat.delimiter("::"));

            return pos;
        }
    }


    RangeRep parseRange(Statement &stat) {
        if(stat.keyword("FULL")) {
            return RangeRep();
        } else {
            PlaceRep first = parsePlace(stat);
            PlaceRep last  = first;
            if(stat.delimiter('/')) last = parsePlace(stat);
            return RangeRep(first, last);
        }
    }


    SRefAttr<double> *parseReference(Statement &stat) {
        std::string objName = parseString(stat, "Object name expected.");

        // Test for attribute name.
        std::string attrName;

        if(stat.delimiter("->")) {
            // Reference to object attribute.
            attrName = parseString(stat, "Attribute name expected.");
        }

        // Parse an index, if present.
        int index = 0;

        if(stat.delimiter('[')) {
            index = int(Round(parseRealConst(stat)));
            parseDelimiter(stat, ']');

            if(index <= 0) {
                throw ParseError("Expressions::parseReference()",
                                 "Index must be positive.");
            }
        }

        // Build reference and fill in the attribute pointer.
        return new SRefAttr<double>(objName, attrName, index);
    }


    TableRowRep parseTableRow(Statement &stat) {
        std::string tabName = parseString(stat, "Table name expected.");
        parseDelimiter(stat, '@');
        return TableRowRep(tabName, parsePlace(stat));
    }


    // ANCILLARY FUNCTIONS.
    // ------------------------------------------------------------------------

    PtrToScalar<double>
    parseTableExpression(Statement &stat, const Table *t) {
        currentTable = t;
        PtrToScalar<double> result = parseReal(stat);
        currentTable = 0;
        return result;
    }


    std::list<Token> parseTokenList(Statement &stat) {
        std::list<Token> result;

        while(!stat.atEnd()) {
            stat.mark();
            Token token = stat.getCurrent();

            // End of list if one of the following tokens occurs outside
            // of brackets.
            if(token.isDel(',') || token.isDel(';') || token.isDel(')') ||
               token.isDel(']') || token.isDel('}')) {
                // Must reread the terminating token.
                stat.restore();
                break;
            }

            result.push_back(token);

            if(token.isDel('(')) {
                parseBracketList(stat, ')', result);
            } else if(token.isDel('[')) {
                parseBracketList(stat, ']', result);
            } else if(token.isDel('{')) {
                parseBracketList(stat, '}', result);
            }
        }

        return result;
    }


    std::vector<std::list<Token> > parseTokenListArray(Statement &stat) {
        std::vector<std::list<Token> > array;

        if(stat.delimiter('{')) {
            // List of token lists within braces.
            do {
                array.push_back(parseTokenList(stat));
            } while(stat.delimiter(','));

            parseDelimiter(stat, '}');
        } else {
            // A single token list.
            array.push_back(parseTokenList(stat));
        }

        return array;
    }


    // LOCAL FUNCTIONS.
    PtrToScalar<bool> parseAnd(Statement &stat) {
        PtrToScalar<bool> result = parseRelation(stat);
        PtrToScalar<bool> second;

        while(stat.delimiter("&&")) {
            second = parseRelation(stat);
            result = SBinary<bool, bool>::make(logAnd, result, second);
        }

        return result;
    }


    PtrToArray<double> parseArrayFactor(Statement &stat) {
        PtrToArray<double> result = parseArrayPrimary(stat);
        PtrToArray<double> second;

        if(stat.delimiter('^')) {
            second = parseArrayPrimary(stat);
            result = new ABinary<double, double>(power, result, second);
        }

        return result;
    }


    PtrToArray<double> parseArrayPrimary(Statement &stat) {
        PtrToArray<double> result;
        double value;

        if(stat.delimiter('(')) {
            result = parseRealArray(stat);
            parseDelimiter(stat, ')');
        } else if(stat.delimiter('{')) {
            ArrayOfPtrs<double> array;
            do {
                array.push_back(parseReal(stat));
            } while(stat.delimiter(','));
            parseDelimiter(stat, '}');
            result = new AList<double>(array);
        } else if(stat.keyword("COLUMN")) {
            result = parseColumnGenerator(stat);
        } else if(stat.keyword("ROW")) {
            result = parseRowGenerator(stat);
        } else if(stat.keyword("TABLE")) {
            result = parseTableGenerator(stat);
        } else if(stat.real(value)) {
            ArrayOfPtrs<double> array;
            array.push_back(new SConstant<double>(value));
            result = new AList<double>(array);
        } else {
            std::string frstName = parseString(stat, "Object name expected.");
            if(stat.delimiter('(')) {
                // Possible function call.
                PtrToArray<double> arg1;
                PtrToArray<double> arg2;
                if(const TFunction2<double, double> *fun = find(table2, frstName)) {
                    arg1 = parseRealArray(stat);
                    parseDelimiter(stat, ',');
                    arg2 = parseRealArray(stat);
                    result = new ABinary<double, double>(*fun, arg1, arg2);
                } else if(const TFunction1<double, double> *fun =
                              find(table1, frstName)) {
                    arg1 = parseRealArray(stat);
                    result = new AUnary<double, double>(*fun, arg1);
                } else if(frstName == "EVAL") {
                    result = parseRealConstArray(stat);
                } else if(const TFunction0<double> *fun = find(table0, frstName)) {
                    ArrayOfPtrs<double> array;
                    array.push_back(SNull<double>::make(*fun));
                    result = new AList<double>(array);
                } else if(const ArrayFun *fun = find(tablea, frstName)) {
                    PtrToArray<double> arg1 = parseRealArray(stat);
                    ArrayOfPtrs<double> array;
                    array.push_back(new ASUnary<double>(*fun, arg1));
                    result = new AList<double>(array);
                } else {
                    throw ParseError("parseArrayPrimary()",
                                     "Invalid array function name \"" + frstName +
                                     "\".");
                }
                parseDelimiter(stat, ')');
            } else if(stat.delimiter("->")) {
                std::string scndName = parseString(stat, "Attribute name expected.");
                result = new ARefExpr<double>(frstName, scndName);
            } else if(stat.delimiter('@')) {
                PlaceRep row = parsePlace(stat);
                // Possible column names.
                if(stat.delimiter("->")) {
                    std::vector<std::string> cols = parseStringArray(stat);
                    result = new ARow(frstName, row, cols);
                } else {
                    throw ParseError("Expressions::parseReal()",
                                     "Expected a column name.");
                }
            } else {
                result = new ARefExpr<double>(frstName, "");
            }
        }

        return result;
    }


    PtrToArray<double> parseArrayTerm(Statement &stat) {
        PtrToArray<double> result = parseArrayFactor(stat);
        PtrToArray<double> second;

        while(true) {
            if(stat.delimiter('*')) {
                second = parseArrayFactor(stat);
                result = new ABinary<double, double>(times, result, second);
            } else if(stat.delimiter('/')) {
                second = parseArrayFactor(stat);
                result = new ABinary<double, double>(divide, result, second);
            } else {
                break;
            }
        }

        return result;
    }


    void parseBracketList(Statement &stat, char close, std::list<Token> &result) {
        while(true) {
            if(stat.atEnd() || stat.delimiter(';')) {
                throw ParseError("Expressions::parseBracketList()",
                                 "Parentheses, brackets or braces do not nest.");
            }

            Token token = stat.getCurrent();
            result.push_back(token);

            if(token.isDel('(')) {
                parseBracketList(stat, ')', result);
            } else if(token.isDel('[')) {
                parseBracketList(stat, ']', result);
            } else if(token.isDel('{')) {
                parseBracketList(stat, '}', result);
            } else if(token.isDel(close)) {
                return;
            }
        }
    }


    PtrToScalar<double> parseFactor(Statement &stat) {
        PtrToScalar<double> result = parsePrimary(stat);
        PtrToScalar<double> second;

        if(stat.delimiter('^')) {
            second = parsePrimary(stat);
            result = SBinary<double, double>::make(power, result, second);
        }

        return result;
    }


    PtrToScalar<double> parsePrimary(Statement &stat) {
        PtrToScalar<double> arg1;
        PtrToScalar<double> arg2;
        PtrToScalar<double> result;
        double value;

        if(stat.delimiter('(')) {
            // Expression in "(...)".
            result = parseReal(stat);
            parseDelimiter(stat, ')');
        } else if(stat.real(value)) {
            // Litteral constant.
            result = new SConstant<double>(value);
        } else if(currentArray.isValid() && stat.delimiter('#')) {
            // "#" hash expression
            result = new SHash(*currentArray);
        } else {
            // Primary beginning with a name.
            std::string frstName = parseString(stat, "Real primary expected.");

            if(stat.delimiter('(')) {
                // Possible function call.
                if(const TFunction2<double, double> *fun = find(table2, frstName)) {
                    arg1 = parseReal(stat);
                    parseDelimiter(stat, ',');
                    arg2 = parseReal(stat);
                    result = SBinary<double, double>::make(*fun, arg1, arg2);
                } else if(const TFunction1<double, double> *fun =
                              find(table1, frstName)) {
                    arg1 = parseReal(stat);
                    result = SUnary<double, double>::make(*fun, arg1);
                } else if(frstName == "EVAL") {
                    result = new SConstant<double>(parseRealConst(stat));
                } else if(const TFunction0<double> *fun = find(table0, frstName)) {
                    result = SNull<double>::make(*fun);
                } else if(const ArrayFun *fun = find(tablea, frstName)) {
                    PtrToArray<double> arg1 = parseRealArray(stat);
                    result = new ASUnary<double>(*fun, arg1);
                } else {
                    throw ParseError("parsePrimary()",
                                     "Unknown function name \"" + frstName + "\".");
                }
                parseDelimiter(stat, ')');
            } else if(stat.delimiter("->")) {
                // Possible parameter or object->attribute clause.
                std::string scndName =
                    parseString(stat, "Attribute or element name expected.");

                // Possible index.
                if(stat.delimiter('[')) {
                    arg2 = parseReal(stat);
                    parseDelimiter(stat, ']');
                    PtrToArray<double> arg1 =
                        new ARefExpr<double>(frstName, scndName);
                    result = new Indexer<double>(arg1, arg2);
                } else {
                    result = new SRefExpr<double>(frstName, scndName);
                }
            } else if(stat.delimiter('@')) {
                PlaceRep row = parsePlace(stat);
                // Possible column name.
                if(stat.delimiter("->")) {
                    std::string col = parseString(stat, "Column name expected.");
                    result = new SCell(frstName, row, col);
                } else {
                    throw ParseError("Expressions::parseReal()",
                                     "Expected a column name.");
                }
            } else {
                // If the name is indexed, the name must be a global vector.
                if(stat.delimiter('[')) {
                    arg2 = parseReal(stat);
                    parseDelimiter(stat, ']');
                    PtrToArray<double> arg1 = new ARefExpr<double>(frstName, "");
                    result = new Indexer<double>(arg1, arg2);
                } else if(currentTable) {
                    // Try making a column expression.
                    result = currentTable->makeColumnExpression(frstName);
                    if(! result.isValid()) {
                        // Could not make a column expression from this name.
                        result = new SRefExpr<double>(frstName, "");
                    }
                } else {
                    // Simple name must be a global variable.
                    result = new SRefExpr<double>(frstName, "");
                }
            }
        }

        return result;
    }


    PtrToScalar<bool> parseRelation(Statement &stat) {
        PtrToScalar<double> first;
        PtrToScalar<double> second;
        PtrToScalar<bool>   result;
        bool value;
        std::string name;

        if(stat.delimiter('(')) {
            result = parseBool(stat);
            parseDelimiter(stat, ')');
        } else if(stat.boolean(value)) {
            result = new SConstant<bool>(value);
        } else {
            stat.mark();
            BoolConstant *bc = 0;
            if(stat.word(name) &&
               (bc = dynamic_cast<BoolConstant *>(OpalData::getInstance()->find(name)))) {
                return new SConstant<bool>(bc->getBool());
            }

            stat.restore();
            first = parseReal(stat);
            const TFunction2<bool, double> *code = 0;
            if(stat.delimiter('<')) {
                code = &less;
            } else if(stat.delimiter('>')) {
                code = &greater;
            } else if(stat.delimiter(">=")) {
                code = &greaterEqual;
            } else if(stat.delimiter("<=")) {
                code = &lessEqual;
            } else if(stat.delimiter("==")) {
                code = &equal;
            } else if(stat.delimiter("!=")) {
                code = &notEqual;
            } else {
                throw ParseError("parseRelation()",
                                 "Invalid boolean expression.");
            }
            second = parseReal(stat);
            result = SBinary<bool, double>::make(*code, first, second);
        }

        return result;
    }


    PtrToScalar<double> parseTerm(Statement &stat) {
        PtrToScalar<double> result = parseFactor(stat);
        PtrToScalar<double> second;

        while(true) {
            if(stat.delimiter('*')) {
                second = parseFactor(stat);
                result = SBinary<double, double>::make(times, result, second);
            } else if(stat.delimiter('/')) {
                second = parseFactor(stat);
                result = SBinary<double, double>::make(divide, result, second);
            } else {
                break;
            }
        }

        return result;
    }


    PtrToArray<double> parseColumnGenerator(Statement &stat) {
        // COLUMN function builds an array from a table column.
        // format: COLUMN(<table>, <column>) or
        //         COLUMN(<table>, <column>, <range>).
        // The word "COLUMN" has already been seen by the caller.
        parseDelimiter(stat, '(');
        std::string tabName = parseString(stat, "Table name expected.");

        // Column name.
        parseDelimiter(stat, ',');
        std::string colName = parseString(stat, "Column name expected.");

        // Optional range specification.
        RangeRep range;
        if(stat.delimiter(',')) {
            range = parseRange(stat);
        }

        parseDelimiter(stat, ')');
        return new AColumn(tabName, colName, range);
    }


    PtrToArray<double> parseRowGenerator(Statement &stat) {
        // ROW function builds an array from a table row.
        // format: ROW(<table>, <place>) or
        //         ROW(<table>, <place>, <columns>).
        // The word "ROW" has already been seen by the caller.
        parseDelimiter(stat, '(');
        std::string tabName = parseString(stat, "Table name expected.");

        // Row position.
        parseDelimiter(stat, ',');
        PlaceRep row = parsePlace(stat);

        // Optional column specifications.
        std::vector<std::string> columns;
        if(stat.delimiter(',')) {
            columns = parseStringArray(stat);
            parseDelimiter(stat, ')');
        }

        return new ARow(tabName, row, columns);
    }

    PtrToArray<double> parseTableGenerator(Statement &stat) {
        // TABLE function generates an array from expressions.
        // format: "TABLE(<n>, <real>)" or
        //         "TABLE(<n1>:<n2>, <real>)" or
        //         "TABLE(<n1>:<n2>:<dn>, <real>)".
        // The word "TABLE" has already been read by the caller.
        parseDelimiter(stat, '(');

        // Read the array index set.
        int frst = int(Round(parseRealConst(stat)));
        int last = 1;
        int step = 1;

        if(stat.delimiter(',')) {
            last = frst;
            frst = 1;
        } else if(stat.delimiter(':')) {
            last = int(Round(parseRealConst(stat)));
            if(stat.delimiter(':')) {
                step = int(Round(parseRealConst(stat)));
            }
            parseDelimiter(stat, ',');
        } else {
            throw ParseError("Expressions::parseTableGenerator()",
                             "Index set incorrect or missing.");
        }

        // Check the array index set.
        if(frst <= 0 || last <= 0 || step <= 0) {
            throw ParseError("Expressions::parseTableGenerator()",
                             "Invalid array index set.");
        }

        // Construct the expression generator.
        currentArray = new ATable(frst, last, step);
        currentArray->defineExpression(parseReal(stat));
        parseDelimiter(stat, ')');

        // The call to release() is required because of type constraints.
        return currentArray.release();
    }

}
