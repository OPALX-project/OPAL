#ifndef OPAL_SNull_HH
#define OPAL_SNull_HH

// ------------------------------------------------------------------------
// $RCSfile: SNull.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SNull<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Expressions/TFunction0.h"
#include "Utilities/DomainError.h"
#include "Utilities/OverflowError.h"
#include <cerrno>
#include <iostream>


namespace Expressions {

    // Class SNull
    // ----------------------------------------------------------------------
    /// A scalar expression without operands.
    //  This expression evaluates a scalar function withoud operands
    //  (e.g. a random generator) and returns the result.

    template <class T>
    class SNull: public Scalar<T> {

    public:

        /// Constructor.
        //  Use an operand-less scalar function.
        explicit SNull(const TFunction0<T> &function);

        SNull(const SNull<T> &);
        virtual ~SNull();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate.
        virtual T evaluate() const;

        /// Make expression.
        //  If the function is not a random generator, evaluate it and
        //  store the result as a constant.
        static Scalar<T> *make(const TFunction0<T> &function);

        /// Print expression.
        virtual void print(std::ostream &str, int precedence = 99) const;

    private:

        // Not implemented.
        SNull();
        void operator=(const SNull &);

        // The operation function.
        const TFunction0<T> &fun;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T> inline
    SNull<T>::SNull(const SNull<T> &rhs):
        Scalar<T>(rhs),
        fun(rhs.fun)
    {}


    template <class T> inline
    SNull<T>::SNull(const TFunction0<T> &function):
        fun(function)
    {}


    template <class T> inline
    SNull<T>::~SNull()
    {}


    template <class T> inline
    Scalar<T> *SNull<T>::clone() const {
        return new SNull<T>(*this);
    }


    template <class T> inline
    T SNull<T>::evaluate() const {
        errno = 0;
        T result = (*fun.function)();

        // Test for run-time evaluation errors.
        switch(errno) {

            case EDOM:
                throw DomainError("SNull:evaluate()");

            case ERANGE:
                // Ignore underflow.
                if(result == T(0)) return result;
                throw OverflowError("SNull:evaluate()");

            default:
                return result;
        }
    }


    template <class T> inline
    Scalar<T> *SNull<T>::make(const TFunction0<T> &fun) {
        return new SNull<T>(fun);
    }


    template <class T> inline
    void SNull<T>::print(std::ostream &stream, int) const {
        stream << fun.name << "()";
    }

}

#endif // OPAL_SNull_HH
