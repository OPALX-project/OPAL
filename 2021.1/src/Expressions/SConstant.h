#ifndef OPAL_SConstant_HH
#define OPAL_SConstant_HH

// ------------------------------------------------------------------------
// $RCSfile: SConstant.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SConstant<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:12:35 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Parser/Token.h"
#include <iomanip>
#include <iostream>
#include <list>


namespace Expressions {

    // Class SConstant
    // ----------------------------------------------------------------------
    /// A scalar constant expression.
    //  This expression returns a scalar constant.

    template <class T> class SConstant: public Scalar<T> {

    public:

        /// Constructor.
        //  Use the value of the constant.
        explicit SConstant(T value);

        virtual ~SConstant();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate.
        virtual T evaluate() const;

        /// Test for constant.
        //  Always true.
        virtual bool isConstant() const;

        /// Print expression.
        virtual void print(std::ostream &str, int precedence) const;

    private:

        // Not implemented.
        SConstant();
        SConstant(const SConstant<T> &);
        void operator=(const SConstant<T> &);

        T value;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T> inline
    SConstant<T>::SConstant(T val):
        value(val)
    {}


    template <class T> inline
    SConstant<T>::~SConstant()
    {}


    template <class T> inline
    Scalar<T> *SConstant<T>::clone() const {
        return new SConstant<T>(value);
    }


    template <class T> inline
    T SConstant<T>::evaluate() const {
        return value;
    }


    template <class T> inline
    bool SConstant<T>::isConstant() const {
        return true;
    }


    // All print() methods must be specialised.
    // ----------------------------------------------------------------------

    template<> inline
    void SConstant<bool>::print(std::ostream &os, int) const {
        os << (value ? "TRUE" : "FALSE");
    }


    template<> inline
    void SConstant<double>::print(std::ostream &os, int) const {
        //std::streamsize old_prec = os.precision(12);
        os << value;
        // ada 15-6-2000 statement unreachable  os.precision(old_prec);
    }


    template<> inline
    void SConstant<std::string>::print(std::ostream &os, int) const {
        os << '"' << value << '"';
    }


    template<> inline
    void SConstant<std::list<Token> >::print(std::ostream &os, int) const {
        for(std::list<Token>::const_iterator token = value.begin();
            token != value.end(); ++token) {
            os << *token;
        }
    }

}

#endif // OPAL_SConstant_HH
