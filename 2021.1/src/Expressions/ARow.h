#ifndef OPAL_ARow_HH
#define OPAL_ARow_HH

// ------------------------------------------------------------------------
// $RCSfile: ARow.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ARow
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/PlaceRep.h"
#include <iosfwd>
#include <vector>


namespace Expressions {

    // Class ARow.
    // ----------------------------------------------------------------------
    /// An array expression defined as a table row.
    //  The result is an array formed by the values found in selected
    //  columns of a given row.

    class ARow: public OArray<double> {

    public:

        /// Constructor.
        //  Identify the table by [b]tabName[/b] and the row by the place
        //  reference [b]row[/b].  The array [b]col[/b] contains the column
        //  names to select the columns by name.
        ARow(const std::string &tabName, const PlaceRep &row,
             const std::vector<std::string> &col);

        ARow(const ARow &);
        ~ARow();

        /// Make clone.
        virtual OArray<double> *clone() const;

        /// Evaluate.
        virtual std::vector<double> evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &os, int precedence = 99) const;

    private:

        // Not implemented.
        ARow();
        const ARow &operator=(const ARow &);

        // The table name.
        const std::string tabName;

        // The row specification.
        PlaceRep position;

        // The std::vector of column names.
        const std::vector<std::string> columns;
    };

}

#endif // OPAL_ARow_HH
