#ifndef OPAL_ABinary_HH
#define OPAL_ABinary_HH

// ------------------------------------------------------------------------
// $RCSfile: ABinary.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: ABinary<T,U>
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/AList.h"
#include "Expressions/TFunction2.h"
#include "Utilities/DomainError.h"
#include "Utilities/OpalException.h"
#include "Utilities/OverflowError.h"
#include <cerrno>
#include <iosfwd>
#include <vector>


namespace Expressions {

    // Template class ABinary
    // ----------------------------------------------------------------------
    /// An array expression with two array operands.
    //  The expression first evaluates the two array operands, and then
    //  loops in parallel over the results and applies the given function
    //  to each pair of components.

    template <class T, class U>
    class ABinary: public AList<T> {

    public:

        /// Constructor.
        //  Use scalar function with two operands and two array expressions.
        ABinary(const TFunction2<T, U> &function,
                PtrToArray<U> left, PtrToArray<U> right);

        ABinary(const ABinary<T, U> &);
        virtual ~ABinary();

        /// Make clone.
        virtual OArray<T> *clone() const;

        /// Evaluate.
        virtual std::vector<T> evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &str, int precedence = 99) const;

    private:

        // Not implemented.
        ABinary();
        void operator=(const ABinary &);

        // The operation object.
        const TFunction2<T, U> &fun;

        // The two operands.
        PtrToArray<U> lft;
        PtrToArray<U> rgt;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T, class U> inline
    ABinary<T, U>::ABinary(const ABinary &rhs):
        AList<T>(), fun(rhs.fun), lft(rhs.lft->clone()), rgt(rhs.rgt->clone())
    {}


    template <class T, class U> inline
    ABinary<T, U>::ABinary(const TFunction2<T, U> &function,
                           PtrToArray<U> left, PtrToArray<U> right):
        AList<T>(), fun(function), lft(left), rgt(right)
    {}


    template <class T, class U> inline
    ABinary<T, U>::~ABinary()
    {}


    template <class T, class U> inline
    OArray<T> *ABinary<T, U>::clone() const {
        return new ABinary<T, U>(*this);
    }


    template <class T, class U> inline
    std::vector<T> ABinary<T, U>::evaluate() const {
        errno = 0;
        std::vector<U> op1 = lft->evaluate();
        std::vector<U> op2 = rgt->evaluate();
        std::vector<T> result;

        if(op1.size() == op2.size()) {
            for(typename std::vector<U>::size_type i = 0; i < op1.size(); ++i) {
                result.push_back((*fun.function)(op1[i], op2[i]));
            }
        } else if(op1.size() == 1) {
            for(typename std::vector<U>::size_type i = 0; i < op2.size(); ++i) {
                result.push_back((*fun.function)(op1[0], op2[i]));
            }
        } else if(op2.size() == 1) {
            for(typename std::vector<U>::size_type i = 0; i < op1.size(); ++i) {
                result.push_back((*fun.function)(op1[i], op2[0]));
            }
        } else {
            throw OpalException("ABinary::evaluate()",
                                "Inconsistent array dimensions.");
        }

        // Test for run-time evaluation errors.
        switch(errno) {

            case EDOM:
                throw DomainError("ABinary::evaluate()");

            case ERANGE:
                throw OverflowError("ABinary::evaluate()");

            default:
                return result;
        }
    }


    template <class T, class U> inline
    void ABinary<T, U>::print(std::ostream &os, int precedence) const {
        if(fun.precedence >= 0) {
            // Binary operation.
            if(fun.precedence <= precedence) os << "(";
            lft->print(os, fun.precedence - 1);
            os << fun.name;
            rgt->print(os, fun.precedence);
            if(fun.precedence <= precedence) os << ")";
        } else {
            // Function.
            os << fun.name << '(';
            lft->print(os, 0);
            os << ',';
            rgt->print(os, 0);
            os << ')';
        }
    }

}

#endif // OPAL_ABinary_HH
