// ------------------------------------------------------------------------
// $RCSfile: List.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: List
//   The class for OPAL LIST commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/List.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Table.h"
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Parser/SimpleStatement.h"
#include "Parser/Token.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <iostream>

// Class List
// ------------------------------------------------------------------------

// The attributes of class List.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        ALL,         // If true, list all columns.
        COLUMN,      // The column names.
        SIZE
    };
}


List::List():
    Action(SIZE, "LIST",
           "The \"LIST\" statement lists a named table.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "LIST");
    itsAttr[ALL] = Attributes::makeBool
                   ("ALL", "Set true to list all columns");
    itsAttr[COLUMN] = Attributes::makeTokenListArray
                      ("COLUMN", "Column specifiers");

    registerOwnership(AttributeHandler::STATEMENT);
}


List::List(const std::string &name, List *parent):
    Action(name, parent)
{}


List::~List()
{}


List *List::clone(const std::string &name) {
    return new List(name, this);
}


void List::execute() {
    std::string tableName = Attributes::getString(itsAttr[TABLE]);
    Table *table = Table::find(tableName);

    if(table) {
        std::string fileName = Attributes::getString(itsAttr[FNAME]);

        if(fileName == "TERM") {
            list(std::cout, table);
        } else {
            std::ofstream os(fileName.c_str());
            if(os.good()) {
                list(os, table);
            } else {
                throw OpalException("List::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw ParseError("List::execute()",
                         "Table \"" + tableName + "\" not found.");
    }
}


void List::list(std::ostream &os, Table *table) {
    bool listAll = Attributes::getBool(itsAttr[ALL]);
    Table::CellArray cells;

    // Column specification desired ?
    if(listAll) {
        cells = table->getDefault();
    } else {
        // Parse the column specifications.
        // These are returned as TokenList's, since the table is not yet known.
        // to the table; it therefore returns the column expressions in the
        // form of a TokenListArray.
        std::vector<std::list<Token> > columns =
            Attributes::getTokenListArray(itsAttr[COLUMN]);

        for(std::vector<std::list<Token> >::iterator i = columns.begin();
            i != columns.end(); ++i) {
            std::list<Token> tokenList = *i;
            SimpleStatement col("COLUMN[]", tokenList);
            col.start();
            int width = 12, prec = 8;
            Expressions::PtrToScalar<double> expr =
                Expressions::parseTableExpression(col, table);

            if(col.delimiter(':')) {
                if(col.integer(width)) {
                    if(col.delimiter(':')) {
                        if(! col.integer(prec)) {
                            throw ParseError("List::list()",
                                             "Invalid <precision> for table column \"" +
                                             getOpalName() + "\".");
                        }
                    }
                } else {
                    throw ParseError("List::list()",
                                     "Invalid <width> for table column \"" +
                                     getOpalName() + "\".");
                }
            }

            if(width < prec + 6) width = prec + 6;
            cells.push_back(Table::Cell(expr, width, prec));
        }
    }

    // Now list the table.
        table->printTable(os, cells);
}