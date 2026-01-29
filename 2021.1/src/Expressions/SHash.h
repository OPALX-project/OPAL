#ifndef OPAL_SHash_HH
#define OPAL_SHash_HH

// ------------------------------------------------------------------------
// $RCSfile: SHash.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SHash
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include <iosfwd>


namespace Expressions {

    // Template class SHash
    // ----------------------------------------------------------------------

    // Forward declaration.
    class ATable;


    /// A scalar expression.
    //  Represents a ``#'' value in a TABLE() expression.

    class SHash: public Scalar<double> {

    public:

        /// Constructor.
        //  Store the TABLE() expression which will return the value of ``#''.
        explicit SHash(const ATable &);

        SHash(const SHash &);
        virtual ~SHash();

        /// Make clone.
        virtual Scalar<double> *clone() const;

        /// Evaluate.
        virtual double evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &str, int precedence = 99) const;

    private:

        // Not implemented.
        SHash();
        void operator=(const SHash &);

        // The ATable object to which this SHash belongs.
        // When  evaluate()  is called, the value for  #  is obtained from
        // "itsTable".
        const ATable &itsTable;
    };

}

#endif // OPAL_SHash_HH
