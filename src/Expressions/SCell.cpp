// ------------------------------------------------------------------------
// $RCSfile: SCell.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class SCell:
//   Reference to a table cell.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/22 15:16:02 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/SCell.h"
#include "AbstractObjects/Table.h"
#include "Channels/ConstChannel.h"

#include <iostream>


// Class SCell
// ------------------------------------------------------------------------


namespace Expressions {

    SCell::SCell
    (const std::string &tab, const PlaceRep &place, const std::string &col):
        tab_name(tab), position(place), col_name(col), itsChannel(0)
    {}


    SCell::~SCell() {
        delete itsChannel;
    }


    Scalar<double> *SCell::clone() const {
        return new SCell(tab_name, position, col_name);
    }


    double SCell::evaluate() const {
        // Define table pointer and refill the table if required.
        Table *table = Table::find(tab_name);
        table->fill();

        // Return value of table cell.
        return table->getCell(position, col_name);
    }


    void SCell::print(std::ostream &os, int) const {
        os << tab_name << "->" << position;
        os << "->" << col_name;
        return;
    }

}
