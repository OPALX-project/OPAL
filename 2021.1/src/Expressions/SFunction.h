#ifndef OPAL_SFunction_HH
#define OPAL_SFunction_HH

// ------------------------------------------------------------------------
// $RCSfile: SFunction.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SFunction
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


// Class SFunction
// ------------------------------------------------------------------------
/// Functions of arc length.
//  This class handles the position funcitions SI(), SC(), and SO().
//  The class SFunction is instantiated by the class AlignHandler, which
//  also takes care of updating the arc length.  The methods arcIn(),
//  arcCtr(), and arcOut() can then return s at the entrance, at the centre,
//  and at the end of the current element respectively.

class SFunction {

public:

    /// Default constructor.
    //  This constructor resets the arc length and registers [b]this[/b]
    //  as the current arc length function.  Only one such function may
    //  be active at any time.
    SFunction();

    /// Destructor.
    //  Unregister [b]this[/b] as the current arc length function.
    ~SFunction();

    /// Return arc length at entrance SI().
    static double arcIn();

    /// Return arc length at center SC().
    static double arcCtr();

    /// Return arc length at exit SO().
    static double arcOut();

    /// Reset the arc length to zero.
    void reset();

    /// Advance position by element length.
    void update(double length);

private:

    // Not implemented.
    SFunction(const SFunction &);
    void operator=(const SFunction &);

    // The function evaluating the arc length.
    double position(double flag) const;

    // The length of the current element.
    double elementLength;

    // The arc length to the exit of the current element.
    double exitArc;

    // The currently active S-function.
    // Only one function may be active at any time.
    static const SFunction *sfun;
};

#endif // OPAL_SFunction_HH
