// ------------------------------------------------------------------------
// $RCSfile: TableRowRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableRowRep
//   A class used to represent a table row reference.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:05:47 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/TableRowRep.h"
#include <iostream>


// Class TableRowRep
// ------------------------------------------------------------------------

TableRowRep::TableRowRep():
    tabName(), position()
{}


TableRowRep::TableRowRep(const TableRowRep &rhs):
    tabName(rhs.tabName), position(rhs.position)
{}


TableRowRep::TableRowRep(const std::string &tab, const PlaceRep &row):
    tabName(tab), position(row)
{}


TableRowRep::~TableRowRep()
{}


const TableRowRep &TableRowRep::operator=(const TableRowRep &rhs) {
    tabName = rhs.tabName;
    position = rhs.position;
    return *this;
}


const std::string &TableRowRep::getTabName() const {
    return tabName;
}


PlaceRep TableRowRep::getPosition() const {
    return position;
}


void TableRowRep::print(std::ostream &os) const {
    os << tabName << "->" << position;
    return;
}
