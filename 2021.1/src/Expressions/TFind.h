#ifndef OPAL_TFind_HH
#define OPAL_TFind_HH

// ------------------------------------------------------------------------
// $RCSfile: TFind.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template function: find<T>(const T table[], const string &name).
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <string>


namespace Expressions {

    // Function find
    // ----------------------------------------------------------------------
    /// Look up name.
    //  The input table is a C array of structures containing a [b]name[/b]
    //  entry.  The result is a pointer to one component of the array or NULL.
    //  May be replaced by a map<> with suitable parameters.

    template <class T> inline
    const T *find(const T table[], const std::string &name) {
        for(const T *x = table; x->name != 0; x++) {
            if(x->name == name) return x;
        }

        return 0;
    }

}

#endif // OPAL_TFind_HH
