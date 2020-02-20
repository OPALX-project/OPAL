#ifndef OPAL_AList_HH
#define OPAL_AList_HH

// ------------------------------------------------------------------------
// $RCSfile: AList.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: AList<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Expressions/SConstant.h"
#include <iostream>
#include <vector>


namespace Expressions {

    // Template Class AList
    // ----------------------------------------------------------------------
    /// An array expression defined by a list of scalar expressions.
    //  Each scalar expression is evaluated separately, and its result is
    //  appended to the result array.

    template <class T>
    class AList: public OArray<T> {

    public:

        /// Default constructor.
        //  Construct empty array.
        AList();

        /// Constructor.
        //  Use an array of scalar expressions.
        explicit AList(const ArrayOfPtrs<T> &);

        AList(const AList &);
        virtual ~AList();

        /// Make clone.
        virtual OArray<T> *clone() const;

        /// Evaluate.
        virtual std::vector<T> evaluate() const;

        /// Print array expression.
        virtual void print(std::ostream &, int precedence = 99) const;

    protected:

        /// The vector of expressions.
        ArrayOfPtrs<T> itsValue;

    private:

        // Not implemented.
        void operator=(const AList &);
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T>
    AList<T>::AList():
        itsValue()
    {}


    template <class T>
    AList<T>::AList(const AList<T> &rhs):
        OArray<T>(rhs),
        itsValue(rhs.itsValue)
    {}


    template <class T>
    AList<T>::AList(const ArrayOfPtrs<T> &value):
        itsValue(value)
    {}


    template <class T>
    AList<T>::~AList()
    {}


    template <class T>
    OArray<T> *AList<T>::clone() const {
        return new AList<T>(*this);
    }


    template <class T>
    std::vector<T> AList<T>::evaluate() const {
        std::vector<T> result(itsValue.size(), T(0));

        for(typename ArrayOfPtrs<T>::size_type i = 0; i < itsValue.size(); ++i) {
            result[i] = itsValue[i]->evaluate();
        }

        return result;
    }


    template <class T>
    void AList<T>::print(std::ostream &os, int /*precedence*/) const {
        os << '{';
        typename ArrayOfPtrs<T>::const_iterator i = itsValue.begin();

        while(i != itsValue.end()) {
            (*i)->print(os, 0);
            if(++i == itsValue.end()) break;
            os << ',';
        }

        os << '}';
    }

}

#endif // OPAL_AList_HH
