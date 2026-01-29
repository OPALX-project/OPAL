#ifndef OPAL_ATable_HH
#define OPAL_ATable_HH

// ------------------------------------------------------------------------
// $RCSfile: ATable.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ATable
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include <iosfwd>
#include <vector>


namespace Expressions {

    // Class ATable.
    // ----------------------------------------------------------------------
    /// An array expression generated from a TABLE() function.
    //  This expression uses one or more AHash objects to represent
    //  the current index in the TABLE expression.  These may retrieve
    //  the index from the ATable object by calling getHash().
    //  The result is computed by setting ``#'' to each index value in
    //  turn, evaluating the scalar expression for this value, and storing
    //  it in the proper place.

    class ATable: public OArray<double> {

    public:

        /// Constructor.
        //  Use the three index arguments from the TABLE() function.
        //  The expression is assigned later by [tt]defineExpression()[/tt].
        ATable(int n1, int n2, int n3);

        ATable(const ATable &);
        ~ATable();

        /// Make clone.
        virtual OArray<double> *clone() const;

        /// Store the generating expression.
        //  A representation of the expression read in the TABLE() function.
        void defineExpression(PtrToScalar<double>);

        /// Evaluate.
        virtual std::vector<double> evaluate() const;

        /// Return the current value of '#'.
        double getHash() const;

        /// Print expression.
        virtual void print(std::ostream &os, int precedence = 99) const;

    private:

        // Not implemented.
        ATable();
        const ATable &operator=(const ATable &);

        // Column name.
        PtrToScalar<double> itsExpr;

        // The range specification.
        int itsBegin, itsEnd, itsStep;

        // The current value of '#'.
        mutable int itsHash;
    };

}

#endif // OPAL_ATable_HH
