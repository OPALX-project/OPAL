// ------------------------------------------------------------------------
// $RCSfile: Table.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Table
//   This abstract base class defines the interface for the OPAL "TABLE"
//   command family.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Table.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Timer.h"
#include <iostream>


// Class Table
// ------------------------------------------------------------------------

Table::~Table() {
    // Remove pointer to "this" from table directory.
    if(! builtin) OpalData::getInstance()->unregisterTable(this);
}


bool Table::canReplaceBy(Object */*newObject*/) {
    return false;
}


Table *Table::find(const std::string &name) {
    Table *table = dynamic_cast<Table *>(OpalData::getInstance()->find(name));
    if(table == 0) {
        throw OpalException("Table::find()", "Table \"" + name + "\" not found.");
    }
    return table;
}


const std::string Table::getCategory() const {
    return "TABLE";
}


bool Table::shouldTrace() const {
    return true;

}

bool Table::shouldUpdate() const {
    return true;
}


void Table::invalidate() {
    if(dynamic) refill = true;
}


Table::Table(int size, const char *name, const char *help):
    Object(size, name, help), dynamic(false), refill(false) {
    // Do not link table exemplar.
}


Table::Table(const std::string &name, Table *parent):
    Object(name, parent), dynamic(true), refill(true) {
    // Link table to table directory.
    OpalData::getInstance()->registerTable(this);
}
