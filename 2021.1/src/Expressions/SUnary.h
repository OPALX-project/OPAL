#ifndef OPAL_SUnary_HH
#define OPAL_SUnary_HH

// ------------------------------------------------------------------------
// $RCSfile: SUnary.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SUnary<T,U>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Expressions/SConstant.h"
#include "Expressions/TFunction1.h"
#include "Utilities/DomainError.h"
#include "Utilities/OverflowError.h"
#include <cerrno>
#include <iostream>


namespace Expressions {

    // Template Class SUnary
    // ----------------------------------------------------------------------
    /// A scalar expression with one scalar operand.
    //  This expression evaluates a scalar operand, and then returns the
    //  result of applying a scalar function to the result.

    template <class T, class U>
    class SUnary: public Scalar<T> {

    public:

        /// Constructor.
        //  Use scalar function of one operand and a scalar operand.
        SUnary(const TFunction1<T, U> &function, PtrToScalar<U> operand);

        SUnary(const SUnary<T, U> &);
        virtual ~SUnary();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate.
        virtual T evaluate() const;

        /// Make a new expression.
        //  If the operand is constant and the function is not a random
        //  generator, evaluate the expression and  store the result as a
        //  constant.
        static Scalar<T> *make(const TFunction1<T, U> &function,
                               PtrToScalar<U> operand);

        /// Print expression.
        virtual void print(std::ostream &, int precedence = 99) const;

    private:

        // Not implemented.
        SUnary();
        void operator=(const SUnary &);

        // The operation object.
        const TFunction1<T, U> &fun;

        // The two operands.
        PtrToScalar<U> opr;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T, class U> inline
    SUnary<T, U>::SUnary(const SUnary<T, U> &rhs):
        Scalar<T>(rhs),
        fun(rhs.fun), opr(rhs.opr->clone())
    {}


    template <class T, class U> inline
    SUnary<T, U>::SUnary(const TFunction1<T, U> &function,
                         PtrToScalar<U> oper):
        fun(function), opr(oper)
    {}


    template <class T, class U> inline
    SUnary<T, U>::~SUnary()
    {}


    template <class T, class U> inline
    Scalar<T> *SUnary<T, U>::clone() const {
        return new SUnary<T, U>(*this);
    }


    template <class T, class U> inline
    T SUnary<T, U>::evaluate() const {
        errno = 0;
        U arg = opr->evaluate();
        T result = (*fun.function)(arg);

        // Test for run-time evaluation errors.
        switch(errno) {

            case EDOM:
                throw DomainError("SUnary::evaluate()");

            case ERANGE:
                // Ignore underflow.
                if(result == T(0)) return result;
                throw OverflowError("SUnary::evaluate()");

            default:
                ;
        }

        return result;
    }


    template <class T, class U> inline
    Scalar<T> *SUnary<T, U>::make
    (const TFunction1<T, U> &function, PtrToScalar<U> operand) {
        // We must pick up the constant flag before the ownership of "operand"
        // is transferred to "result".
        bool isConst = operand->isConstant();
        PtrToScalar<T> result = new SUnary<T, U>(function, operand);;

        if(function.precedence != -2) {
            if(isConst) {
                // Constant expression.
                double value = result->evaluate();
                result = new SConstant<T>(value);
            }
        }

        // Other expression.
        return result.release();
    }


    template <class T, class U> inline
    void SUnary<T, U>::print(std::ostream &os, int precedence) const {
        if(fun.precedence >= 0) {
            if(fun.precedence <= precedence) os << "(";
            os << fun.name;
            opr->print(os, fun.precedence);
            if(fun.precedence <= precedence) os << ")";
        } else {
            os << fun.name << "(";
            opr->print(os, 0);
            os << ")";
        }
    }

}

#endif // OPAL_SUnary_HH
