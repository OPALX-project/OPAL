#ifndef OPAL_SValue_HH
#define OPAL_SValue_HH

// ------------------------------------------------------------------------
// $RCSfile: SValue.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SValue<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#include "Parser/Token.h"
#include <iosfwd>
#include <iostream>
#include <list>


namespace Expressions {

    // Template class SValue
    // ----------------------------------------------------------------------
    /// Object attribute with a constant scalar value.

    template <class T> class SValue: public AttributeBase {

    public:

        /// Default constructor.
        //  Construct zero value.
        SValue();

        /// Constructor.
        //  Use constant value.
        explicit SValue(const T &val);

        SValue(const SValue<T> &);
        virtual ~SValue();

        /// Make clone.
        virtual SValue<T> *clone() const;

        /// Evaluate.
        //  Return the (already known) value.
        virtual T evaluate();

        /// Print the attribute value.
        virtual void print(std::ostream &) const;

    protected:

        /// The value of the attribute.
        mutable T value;

    private:

        // Not implemented.
        void operator=(const SValue<T> &);
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T>
    SValue<T>::SValue():
        value(T(0))
    {}


    template <class T>
    SValue<T>::SValue(const SValue<T> &rhs):
        value(rhs.value)
    {}


    template <class T>
    SValue<T>::SValue(const T &val):
        value(val)
    {}


    template <class T>
    SValue<T>::~SValue()
    {}


    template <class T>
    SValue<T> *SValue<T>::clone() const {
        return new SValue<T>(value);
    }


    template <class T>
    T SValue<T>::evaluate() {
        return value;
    }


    // Print methods are specialised.
    // ------------------------------------------------------------------------

    template <class T>
    void SValue<T>::print(std::ostream &os) const {
        os << value;
        return;
    }


    template<> inline
    void SValue<std::list<Token> >::print(std::ostream &os) const {
        for(std::list<Token>::iterator token = value.begin();
            token != value.end(); ++token) {
            os << *token;
        }

        return;
    }


    template<> inline
    void SValue<bool>::print(std::ostream &os) const {
        os << (value ? "TRUE" : "FALSE");
        return;
    }


    template<> inline
    void SValue<std::string>::print(std::ostream &os) const {
        os << '"' << value << '"';
        return;
    }

}

#endif // OPAL_SValue_HH
