// ------------------------------------------------------------------------
// $RCSfile: ARow.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ARow
//   Expression class. Generates a list of values from a TABLE() function.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/22 15:16:02 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/ARow.h"
#include "AbstractObjects/Table.h"

#include <iostream>


// Class ARow
// ------------------------------------------------------------------------

namespace Expressions {

    ARow::ARow(const ARow &rhs):
        OArray<double>(rhs), tabName(rhs.tabName), position(rhs.position),
        columns(rhs.columns)
    {}


    ARow::ARow(const std::string &tName, const PlaceRep &row,
               const std::vector<std::string> &cols):
        OArray<double>(), tabName(tName), position(row), columns(cols)
    {}


    ARow::~ARow()
    {}


    OArray<double> *ARow::clone() const {
        return new ARow(*this);
    }


    std::vector<double> ARow::evaluate() const {
        Table *table = Table::find(tabName);
        table->fill();
        return table->getRow(position, columns);
    }


    void ARow::print(std::ostream &os, int) const {
        os << "ROW(" << tabName << ',' << position << ",{";
        std::vector<std::string>::const_iterator i = columns.begin();

        while(i != columns.end()) {
            os << '"' << *i << '"';
            if(++i == columns.end()) break;
            os << ',';
        }
        os << "})";
    }

}
