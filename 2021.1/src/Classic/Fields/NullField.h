#ifndef CLASSIC_NullField_HH
#define CLASSIC_NullField_HH

// ------------------------------------------------------------------------
// $RCSfile: NullField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NullField
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


// Class NullField
// ------------------------------------------------------------------------
/// A zero electromagnetic field.

class NullField: public EMField {

public:

    NullField();
    NullField(const NullField &right);
    virtual ~NullField();
    const NullField &operator=(const NullField &right);

    /// Scale the field.
    //  Multiply the field by [b]scalar[/b].
    //  Obviously this method does nothing.
    virtual void scale(double scalar);
};

#endif // CLASSIC_NullField_HH
