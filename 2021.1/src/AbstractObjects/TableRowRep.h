#ifndef OPAL_TableRowRep_HH
#define OPAL_TableRowRep_HH

// ------------------------------------------------------------------------
// $RCSfile: TableRowRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableRowRep
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/PlaceRep.h"
#include <string>

class Table;


// Class TableRowRep
// ------------------------------------------------------------------------
/// Representation of a table row reference.
//  Such a reference consists of two parts:
//  [ol]
//  [li]A table name,
//  [li]A place reference which selects the table row.
//  [/ol]

class TableRowRep {

public:

    /// Default constructor.
    //  Constructs undefined reference.
    TableRowRep();

    /// Constructor.
    //  Construct reference to the row identified by [b]row[/b] of the
    //  table with name [b]name[/b].
    TableRowRep(const std::string &tab, const PlaceRep &row);

    TableRowRep(const TableRowRep &);
    ~TableRowRep();
    const TableRowRep &operator=(const TableRowRep &);

    /// Return the table name.
    const std::string &getTabName() const;

    /// Return the row position representation.
    PlaceRep getPosition() const;

    /// Print in input format.
    void print(std::ostream &os) const;

private:

    // Table row structure.
    std::string tabName;      // Table name.
    PlaceRep position;   // Position reference.
};


inline std::ostream &operator<<(std::ostream &os, const TableRowRep &row) {
    row.print(os);
    return os;
}

#endif // OPAL_TableRowRep_HH
