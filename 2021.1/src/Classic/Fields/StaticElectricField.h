#ifndef CLASSIC_StaticElectricField_HH
#define CLASSIC_StaticElectricField_HH

// ------------------------------------------------------------------------
// $RCSfile: StaticElectricField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StaticElectricField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:36 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/EMField.h"


// Class StaticElectricField
// ------------------------------------------------------------------------
/// Abstract base class for static electric fields.

class StaticElectricField: public EMField {

public:

    StaticElectricField();
    virtual ~StaticElectricField();
};

#endif // CLASSIC_StaticElectricField_HH
