#ifndef OPAL_ASUnary_HH
#define OPAL_ASUnary_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ASUnary.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: ASUnary<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Expressions/SConstant.h"
#include "Expressions/TFunction1.h"
#include "Utilities/DomainError.h"
#include "Utilities/OverflowError.h"
#include <cerrno>
#include <iosfwd>
#include <vector>


namespace Expressions {

    // Template Class ASUnary
    // ----------------------------------------------------------------------
    /// A scalar expression with one array operand.
    //  This expression first evaluates its operand to get an array.
    //  The scalar result is computed then by applying a scalar-valued
    //  function to this array (e. g. minimum or maximum).

    template <class T>
    class ASUnary: public Scalar<T> {

    public:

        /// Constructor.
        //  Use a function which maps a std::vector<T> to a T and an array
        //  expression.
        ASUnary(const TFunction1<T, const std::vector<T>&> &function,
                PtrToArray<T> &oper);

        ASUnary(const ASUnary<T> &);
        virtual ~ASUnary();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate expression.
        virtual T evaluate() const;

        /// Make new expression.
        //  If all component of the argument are constant, then evaluate
        //  and store the result as a constant.
        static PtrToScalar<T> make
        (const TFunction1<T, const std::vector<T>&> &function, PtrToArray<T> &oper);

        /// Print expression.
        virtual void print(std::ostream &, int precedence = 99) const;

    private:

        // Not implemented.
        ASUnary();
        void operator=(const ASUnary &);

        // The operation object.
        const TFunction1<T, const std::vector<T>&> &fun;

        // The operand.
        PtrToArray<T> opr;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T> inline
    ASUnary<T>::ASUnary(const ASUnary<T> &rhs):
        Scalar<T>(rhs),
        fun(rhs.fun), opr(rhs.opr->clone())
    {}


    template <class T> inline
    ASUnary<T>::ASUnary(const TFunction1<T, const std::vector<T>&> &function,
                        PtrToArray<T> &oper):
        fun(function), opr(oper)
    {}


    template <class T> inline
    ASUnary<T>::~ASUnary()
    {}


    template <class T> inline
    Scalar<T> *ASUnary<T>::clone() const {
        return new ASUnary<T>(*this);
    }


    template <class T> inline
    T ASUnary<T>::evaluate() const {
        errno = 0;
        const std::vector<T> arg = opr->evaluate();
        T result = (*fun.function)(arg);

        // Test for run-time evaluation errors.
        switch(errno) {

            case EDOM:
                throw DomainError("ASUnary::evaluate()");

            case ERANGE:
                // Ignore underflow.
                if(result == T(0)) return result;
                throw OverflowError("ASUnary::evaluate()");

            default:
                ;
        }

        return result;
    }


    template <class T> inline
    PtrToScalar<T> ASUnary<T>::make
    (const TFunction1<T, const std::vector<T>&> &function, PtrToArray<T> &oper) {
        // We must pick up the constant flag before the ownership of "oper"
        // is transferred to "result".
        bool isConst = oper->isConstant();
        PtrToScalar<T> result = new ASUnary<T>(function, oper);

        if(function.precedence != -2) {
            // Not a random expression.
            if(isConst) {
                // Replace constant expression by its value.
                result = new SConstant<T>(result->evaluate());
            }
        }

        // Other expression.
        return result.release();
    }


    template <class T> inline
    void ASUnary<T>::print(std::ostream &os, int precedence) const {
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

#endif // OPAL_ASUnary_HH
