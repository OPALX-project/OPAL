#ifndef CLASSIC_StaticMagneticField_HH
#define CLASSIC_StaticMagneticField_HH

// ------------------------------------------------------------------------
// $RCSfile: StaticMagneticField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StaticMagneticField
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


// Class StaticMagneticField
// ------------------------------------------------------------------------
/// Abstract base class for static magnetic fields.

class StaticMagneticField: public EMField {

public:

    StaticMagneticField();
    virtual ~StaticMagneticField();
};

#endif // CLASSIC_StaticMagneticField_HH
