// ------------------------------------------------------------------------
// $RCSfile: TableRow.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableRow
//   A class used to parse a reference to a table lines.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/TableRow.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using namespace Expressions;


// Class TableRow
// ------------------------------------------------------------------------

namespace Attributes {

    TableRow::TableRow(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    TableRow::~TableRow()
    {}


    const std::string &TableRow::getType() const {
        static const std::string type("table line");
        return type;
    }


    void TableRow::parse(Attribute &attr, Statement &stat, bool) const {
        attr.set(new SValue<TableRowRep>(parseTableRow(stat)));
    }

};
