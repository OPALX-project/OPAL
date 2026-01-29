#ifndef OPAL_TFunction2_HH
#define OPAL_TFunction2_HH

// ------------------------------------------------------------------------
// $RCSfile: TFunction2.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: TFunction2<T,U>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


namespace Expressions {

    // Class TFunction2
    // ----------------------------------------------------------------------
    /// A function of two U's returning a T.
    //  This structure groups the properties of a function with two arguments,
    //  namely its name, its precedence, and a pointer to the function to
    //  evaluate it.

    template <class T, class U> struct TFunction2 {

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
        T(*function)(U, U);
    };

}

#endif // OPAL_TFunction2_HH
