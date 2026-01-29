#ifndef OPAL_List_HH
#define OPAL_List_HH

// ------------------------------------------------------------------------
// $RCSfile: List.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: List
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "AbstractObjects/Attribute.h"
#include <iosfwd>

class Table;


// Class List
// ------------------------------------------------------------------------
/// The LIST command.

class List: public Action {

public:

    /// Exemplar constructor.
    List();

    virtual ~List();

    /// Make clone.
    virtual List *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    List(const List &);
    void operator=(const List &);

    // Clone constructor.
    List(const std::string &name, List *parent);

    // Do the listing.
    void list(std::ostream &, Table *);
};

#endif // OPAL_List_HH
