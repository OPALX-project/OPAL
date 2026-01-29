#ifndef OPAL_SCell_HH
#define OPAL_SCell_HH

// ------------------------------------------------------------------------
// $RCSfile: SCell.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class SCell:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/PlaceRep.h"
#include <iosfwd>

class ConstChannel;


namespace Expressions {

    // Class SCell
    // ----------------------------------------------------------------------
    /// A scalar expression referring to a table cell.
    //  This expression returns the value in a table cell, identified by
    //  the table name, the row (a place reference), and the column (a name).

    class SCell: public Scalar<double> {

    public:

        /// Constructor.
        //  Identify the table by its name [b]tab[/b], the row by the place
        //  reference [b]place[/b] and the column by its name [b]col[/b].
        SCell(const std::string &tab, const PlaceRep &place, const std::string &col);

        virtual ~SCell();

        /// Make clone.
        virtual Scalar<double> *clone() const;

        /// Evaluate.
        //  Evaluate the reference and return the value in that cell.
        virtual double evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &stream, int) const;

    private:

        // Not implemented.
        SCell();
        void operator=(const SCell &);

        // Names of table, of position, and of column, plus occurrence count.
        const std::string tab_name;
        const PlaceRep position;
        const std::string col_name;

        // The Chanel leading to the table cell.
        mutable ConstChannel *itsChannel;
    };

}

#endif // OPAL_SCell_HH
