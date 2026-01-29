#ifndef OPAL_Table_HH
#define OPAL_Table_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Table.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Table
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Object.h"
#include "AbstractObjects/Expressions.h"
#include <vector>

class Attribute;
class Beamline;
class PlaceRep;
class RangeRep;


// Class Table
// ------------------------------------------------------------------------
/// The base class for all OPAL tables.
//  It implements the common behaviour of tables, it can also be used via
//  dynamic casting to determine whether an object represents a table.
//  During matching operations, a Table must often be recalculated.  This
//  is done by calling its [tt]fill()[/tt] method.  To avoid unnecessary
//  filling, the user may define a table as ``static'' in the command
//  creating the table.  In this case OPAL assumes that the table never
//  changes.

class Table: public Object {

public:

    /// Descriptor for printing a table cell.
    struct Cell {
        /// Constructor.
        Cell(Expressions::PtrToScalar<double> expr, int width, int prec):
            itsExpr(expr), printWidth(width), printPrecision(prec) {}

        /// The expression generating the values for this Cell.
        Expressions::PtrToScalar<double> itsExpr;

        /// The cell width in print characters.
        int printWidth;

        /// The number of digits printed for this cell.
        int printPrecision;
    };

    /// An array of cell descriptors.
    typedef std::vector<Cell> CellArray;

    virtual ~Table();

    /// Test if object can be replaced.
    //  A table cannot be replaced, thus this method never returns true.
    virtual bool canReplaceBy(Object *newObject);

    /// Refill the buffer.
    //  If [b]refill[/b] is true, call the table algorithm to fill the buffer.
    virtual void fill() = 0;

    /// Find named Table.
    //  If a table with name [b]name[/b] exists, return a pointer to that table.
    //  If no such table exists, throw OpalException.
    static Table *find(const std::string &name);

    /// Return the object category as a string.
    //  Return the string "TABLE".
    virtual const std::string getCategory() const;

    /// Trace flag.
    //  If true, the object's execute() function should be traced.
    //  Always true for tables.
    virtual bool shouldTrace() const;

    /// Update flag.
    //  If true, the data structure should be updated before calling execute().
    //  Always true for tables.
    virtual bool shouldUpdate() const;


    /// Return value in selected table cell.
    //  Return the value stored in the table cell identified by the row
    //  at [b]row[/b] and the column name [b]col[/b].
    virtual double getCell(const PlaceRep &row, const std::string &col) = 0;

    /// Return column [b]col[/b] of this table, limited by [b]range[/b].
    virtual std::vector<double> getColumn
    (const RangeRep &range, const std::string &col) = 0;

    /// Return the default print columns.
    //  Returns an array of column descriptors, which, when applied to a
    //  row of the table, gives the default print columns for this table.
    virtual CellArray getDefault() const = 0;

    /// Return the length of the table.
    //  Returns the geometric length of the underlying beam line.
    virtual double getLength() = 0;

    /// Return embedded CLASSIC beamline.
    //  Returns the CLASSIC beamline representing the table.
    //  The data of the table are attached to each position in the line.
    virtual const Beamline *getLine() const = 0;

    /// Return a table row.
    //  Returns the values stored in the row specified by the first argument,
    //  with columns selected by the names stored in the second argument.
    virtual std::vector<double>
    getRow(const PlaceRep &, const std::vector<std::string> &) = 0;

    /// Mark this table as invalid, if it is dynamic.
    virtual void invalidate();

    /// Find out if table depends on the object identified by [b]name[/b].
    //  Must be overridden in derived classes.
    virtual bool isDependent(const std::string &name) const = 0;

    virtual Expressions::PtrToScalar<double>
    makeColumnExpression(const std::string &) const = 0;

    /// Check that [b]rhs[/b] is of same type as [b]this[/b].
    virtual bool matches(Table *rhs) const = 0;

    /// Print list for the table.
    virtual void printTable(std::ostream &, const CellArray &) const = 0;


protected:

    /// Constructor for exemplars.
    Table(int size, const char *name, const char *help);

    /// Constructor for clones.
    Table(const std::string &name, Table *parent);


    /// Flag dynamic table.
    //  If true, the table is dynamic.  If it is also invalidated,
    //  it will be refilled when fill() is called.
    bool dynamic;

    /// Refill flag.
    //  If true, the table has been invalidated.  If it is also dynammic,
    //  it will be refilled when fill() is called.
    bool refill;

private:

    // Not implemented.
    Table();
    Table(const Table &);
    void operator=(const Table &);
};

#endif // OPAL_Table_HH
