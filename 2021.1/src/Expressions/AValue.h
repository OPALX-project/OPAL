#ifndef OPAL_AValue_HH
#define OPAL_AValue_HH

// ------------------------------------------------------------------------
// $RCSfile: AValue.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: AValue<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#include "Expressions/SConstant.h"
#include "Parser/Token.h"
#include <iosfwd>
#include <iostream>
#include <vector>


namespace Expressions {

    // Template class AValue
    // ----------------------------------------------------------------------
    /// Object attribute with a constant array value.

    template <class T> class AValue: public AttributeBase {

    public:

        /// Default constructor.
        //  Construct empty array.
        AValue();

        /// Constructor.
        //  Use the vector of values.
        explicit AValue(const std::vector<T> &val);

        AValue(const AValue &);
        virtual ~AValue();

        /// Make clone.
        virtual AValue *clone() const;

        /// Evaluate.
        //  Return the (already known) value.
        virtual std::vector<T> evaluate();

        /// Print the attribute value.
        virtual void print(std::ostream &) const;

    protected:

        /// The value of the attribute.
        mutable std::vector<T> value;

    private:

        // Not implemented.
        void operator=(const AValue<T> &);
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T>
    AValue<T>::AValue():
        value()
    {}


    template <class T>
    AValue<T>::AValue(const AValue<T> &rhs):
        value(rhs.value)
    {}


    template <class T>
    AValue<T>::AValue(const std::vector<T> &val):
        value(val)
    {}


    template <class T>
    AValue<T>::~AValue()
    {}


    template <class T>
    AValue<T> *AValue<T>::clone() const {
        return new AValue<T>(value);
    }


    template <class T>
    std::vector<T> AValue<T>::evaluate() {
        return value;
    }


    // Print methods are specialised.
    // ------------------------------------------------------------------------

    template<> inline
    void AValue<bool>::print(std::ostream &os) const {
        os << '{';
        std::vector<bool>::const_iterator i = value.begin();

        while(i != value.end()) {
            os << (*i ? "TRUE" : "FALSE");
            if(++i == value.end()) break;
            os << ',';
        }

        os << '}';
        return;
    }


    template<> inline
    void AValue<double>::print(std::ostream &os) const {
        std::streamsize old_prec = os.precision(12);
        os << '{';
        std::vector<double>::const_iterator i = value.begin();

        while(i != value.end()) {
            os << *i;
            if(++i == value.end()) break;
            os << ',';
        }

        os << '}';
        os.precision(old_prec);
        return;
    }


    template<> inline
    void AValue<std::string>::print(std::ostream &os) const {
        os << '{';
        std::vector<std::string>::const_iterator i = value.begin();

        while(i != value.end()) {
            os << '"' << *i << '"';
            if(++i == value.end()) break;
            os << ',';
        }

        os << '}';
        return;
    }


    template<> inline
    void AValue<std::list<Token> >::print(std::ostream &os) const {
        os << '{';
        std::vector<std::list<Token> >::const_iterator i = value.begin();

        while(i != value.end()) {
            for(std::list<Token>::const_iterator j = i->begin();
                j != i->end(); ++j) {
                os << *j;
            }

            if(++i == value.end()) break;
            os << ',';
        }

        os << '}';
        return;
    }

}

#endif // OPAL_AValue_HH
