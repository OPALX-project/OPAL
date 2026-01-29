#ifndef OPAL_AUnary_HH
#define OPAL_AUnary_HH

// ------------------------------------------------------------------------
// $RCSfile: AUnary.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: AUnary<T,U>
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/AList.h"
#include "Expressions/TFunction1.h"
#include "Utilities/DomainError.h"
#include "Utilities/OverflowError.h"
#include <cerrno>
#include <iostream>
#include <vector>


namespace Expressions {

    // Template Class AUnary
    // ----------------------------------------------------------------------
    /// An array expression with one array operand.
    //  This expression first computes its array operand, and then applies
    //  the scalar function to each component and stores the results in the
    //  result vector.

    template <class T, class U>
    class AUnary: public AList<T> {

    public:

        /// Constructor.
        //  Use scalar function of one operand and an array.
        AUnary(const TFunction1<T, U> &function, PtrToArray<U> operand);

        AUnary(const AUnary<T, U> &);
        virtual ~AUnary();

        /// Make clone.
        virtual OArray<T> *clone() const;

        /// Evaluate.
        virtual std::vector<T> evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &, int precedence = 99) const;

    private:

        // Not implemented.
        AUnary();
        void operator=(const AUnary &);

        // The operation object.
        const TFunction1<T, U> &fun;

        // The operand.
        PtrToArray<U> opr;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T, class U> inline
    AUnary<T, U>::AUnary(const AUnary<T, U> &rhs):
        AList<T>(), fun(rhs.fun), opr(rhs.opr->clone())
    {}


    template <class T, class U> inline
    AUnary<T, U>::AUnary(const TFunction1<T, U> &function, PtrToArray<U> oper):
        AList<T>(), fun(function), opr(oper)
    {}


    template <class T, class U> inline
    AUnary<T, U>::~AUnary()
    {}


    template <class T, class U> inline
    OArray<T> *AUnary<T, U>::clone() const {
        return new AUnary<T, U>(*this);
    }


    template <class T, class U> inline
    std::vector<T> AUnary<T, U>::evaluate() const {
        errno = 0;
        std::vector<U> arg = opr->evaluate();
        std::vector<T> result;

        for(typename std::vector<U>::size_type i = 0; i < arg.size(); ++i) {
            result.push_back((*fun.function)(arg[i]));
        }

        // Test for run-time evaluation errors.
        switch(errno) {

            case EDOM:
                throw DomainError("AUnary::evaluate()");

            case ERANGE:
                throw OverflowError("AUnary::evaluate()");

            default:
                ;
        }

        return result;
    }


    template <class T, class U> inline
    void AUnary<T, U>::print(std::ostream &os, int precedence) const {
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

#endif // OPAL_AUnary_HH
