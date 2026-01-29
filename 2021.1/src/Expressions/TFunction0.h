#ifndef OPAL_TFunction0_HH
#define OPAL_TFunction0_HH

// ------------------------------------------------------------------------
// $RCSfile: TFunction0.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: TFunction0
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


namespace Expressions {

    // Class TFunction0
    // ----------------------------------------------------------------------
    /// An operand-less function returning a T.
    //  This structure groups the properties of a function without arguments
    //  (e.g. a random generator), namely its name, its precedence, and a
    //  pointer to the function to evaluate it.

    template <class T> struct TFunction0 {

        /// The function name or operator representation.
        const char *name;

        /// The operator precedence.
        //  May be one of the following:
        //  [dl]
        //  [dt]-2[dd] random generator,
        //  [dt]-1[dd] ordinary function,
        //  [dt]other[dd] operators.
        int precedence;

        /// The actual operation.
        //  A pointer to function returning a T.
        T(*function)();
    };

}

#endif // OPAL_TFunction0_HH
