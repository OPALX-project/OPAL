#ifndef OPAL_ObjectFunction_HH
#define OPAL_ObjectFunction_HH

// ------------------------------------------------------------------------
// $RCSfile: ObjectFunction.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ObjectFunction
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

class Object;


// Class ObjectFunction
// ------------------------------------------------------------------------
/// Abstract base class for functor objects whose argument is an Object.
//  The classes derived from this object are normally used to apply an
//  operation to all objects in the OPAL directory.  This is done by calling
//  Directory::apply(functor).  The functor itself takes a pointer to an
//  object as an argument, and defines an operator with the signature
//  [tt]void operator()(Object *)[/tt] which performs that operation.

struct ObjectFunction {

    virtual ~ObjectFunction();

    /// The function to be executed.
    virtual void operator()(Object *) const = 0;
};

#endif // OPAL_ObjectFunction_HH
