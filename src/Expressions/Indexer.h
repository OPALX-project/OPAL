#ifndef OPAL_Indexer_HH
#define OPAL_Indexer_HH

// ------------------------------------------------------------------------
// $RCSfile: Indexer.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Indexer<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "Utilities/CLRangeError.h"
#include "Utilities/OpalException.h"
#include "Utilities/Round.h"
#include <cerrno>
#include <iosfwd>
#include <vector>


namespace Expressions {

    // Template class Indexer
    // ----------------------------------------------------------------------
    /// A scalar expression to retrieve an indexed component from an array.
    //  This expression first evaluates the array operand expression and the
    //  scalar index expression.  Then it applies the index to the array and
    //  returns the selected component.

    template <class T>
    class Indexer: public Scalar<T> {

    public:

        /// Constructor.
        //  Use array expression and index into the array.
        Indexer(PtrToArray<T> array, PtrToScalar<double> index);

        Indexer(const Indexer<T> &);
        virtual ~Indexer();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate.
        virtual T evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &str, int precedence = 99) const;

    private:

        // Not implemented.
        Indexer();
        void operator=(const Indexer &);

        // The two operands.
        PtrToArray<T> lft;
        PtrToScalar<double> rgt;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T> inline
    Indexer<T>::Indexer(const Indexer &rhs):
        Scalar<T>(rhs),
        lft(rhs.lft->clone()), rgt(rhs.rgt->clone())
    {}


    template <class T> inline
    Indexer<T>::Indexer(PtrToArray<T> left, PtrToScalar<double> right):
        lft(left), rgt(right)
    {}


    template <class T> inline
    Indexer<T>::~Indexer()
    {}


    template <class T> inline
    Scalar<T> *Indexer<T>::clone() const {
        return new Indexer<T>(*this);
    }


    template <class T> inline
    T Indexer<T>::evaluate() const {
        std::vector<T> op1 = lft->evaluate();
        int op2 = int(Round(rgt->evaluate()));

        if(op2 <= 0 || static_cast<unsigned>(op2) > op1.size()) {
            throw CLRangeError("Expressions::Indexer()", "Index out of range.");
        }

        return op1[op2-1];
    }


    template <class T> inline
    void Indexer<T>::print(std::ostream &os, int /*precedence*/) const {
        lft->print(os, 0);
        os << '[';
        rgt->print(os, 0);
        os << ']';
    }

}

#endif // OPAL_Indexer_HH
