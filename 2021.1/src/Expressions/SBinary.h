#ifndef OPAL_SBinary_HH
#define OPAL_SBinary_HH

// ------------------------------------------------------------------------
// $RCSfile: SBinary.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SBinary<T,U>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Expressions/SConstant.h"
#include "Expressions/TFunction2.h"
#include "Utilities/DomainError.h"
#include "Utilities/OverflowError.h"
#include <cerrno>
#include <iosfwd>


namespace Expressions {

    // Template class SBinary
    // ----------------------------------------------------------------------
    /// A scalar expression with two scalar operands.
    //  The expression first evaluates the two scalar operands, and then
    //  it applies the given function to the operands and returns the result.

    template <class T, class U>
    class SBinary: public Scalar<T> {

    public:

        /// Constructor.
        //  Use scalar function with two operands and two scalars.
        SBinary(const TFunction2<T, U> &function,
                PtrToScalar<U> left, PtrToScalar<U> right);

        SBinary(const SBinary<T, U> &);
        virtual ~SBinary();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate.
        virtual T evaluate() const;

        /// Make a new expression.
        //  If both operands are constants and the function is not a random
        //  generator, evaluate the expression and  store the result as a
        //  constant.
        static Scalar<T> *make(const TFunction2<T, U> &,
                               PtrToScalar<U> left, PtrToScalar<U> right);

        /// Print expression.
        virtual void print(std::ostream &str, int precedence = 99) const;

    private:

        // Not implemented.
        SBinary();
        void operator=(const SBinary &);

        // The operation object.
        const TFunction2<T, U> &fun;

        // The two operands.
        PtrToScalar<U> lft;
        PtrToScalar<U> rgt;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T, class U> inline
    SBinary<T, U>::SBinary(const SBinary &rhs):
        Scalar<T>(rhs),
        fun(rhs.fun), lft(rhs.lft->clone()), rgt(rhs.rgt->clone())
    {}


    template <class T, class U> inline
    SBinary<T, U>::SBinary(const TFunction2<T, U> &function,
                           PtrToScalar<U> left, PtrToScalar<U> right):
        fun(function), lft(left), rgt(right)
    {}


    template <class T, class U> inline
    SBinary<T, U>::~SBinary()
    {}


    template <class T, class U> inline
    Scalar<T> *SBinary<T, U>::clone() const {
        return new SBinary<T, U>(*this);
    }


    template <class T, class U> inline
    T SBinary<T, U>::evaluate() const {
        errno = 0;
        U op1 = lft->evaluate();
        U op2 = rgt->evaluate();
        T result = (*fun.function)(op1, op2);

        // Test for run-time evaluation errors.
        switch(errno) {

            case EDOM:
                throw DomainError("SBinary::evaluate()");

            case ERANGE:
                // Ignore underflow.
                if(result == T(0)) return result;
                throw OverflowError("SBinary::evaluate()");

            default:
                return result;
        }
    }


    template <class T, class U> inline
    Scalar<T> *SBinary<T, U>::make
    (const TFunction2<T, U> &function,
     PtrToScalar<U> left, PtrToScalar<U> right) {
        // We must pick up the constant flag before the ownerships of "left" and
        // "right" are transferred to "result".
        bool isConst = left->isConstant() && right->isConstant();
        PtrToScalar<T> result = new SBinary<T, U>(function, left, right);

        if(function.precedence != -2) {
            if(isConst) {
                // Replace constant expression by its value.
                result = new SConstant<T>(result->evaluate());
            }
        }

        // Other expression.
        return result.release();
    }


    template <class T, class U> inline
    void SBinary<T, U>::print(std::ostream &os, int precedence) const {
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

        return;
    }

}

#endif // OPAL_SBinary_HH
