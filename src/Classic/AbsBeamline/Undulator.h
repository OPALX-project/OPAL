#ifndef CLASSIC_Undulator_HH
#define CLASSIC_Undulator_HH

// ------------------------------------------------------------------------
// $RCSfile: Undulator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Undulator
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"



// Class Undulator
// ------------------------------------------------------------------------
/// Interface for drift space.
//  Class Undulator defines the abstract interface for a drift space.

class Undulator: public Component {

public:

    /// Constructor with given name.
    explicit Undulator(const std::string &name);

    Undulator();
    Undulator(const Undulator &right);
    virtual ~Undulator();

    /// Apply visitor to Undulator.
    virtual void accept(BeamlineVisitor &) const;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void finalise();

    virtual bool bends() const;

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    //set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices); // Philippe was here

    //set number of slices for map tracking
    std::size_t getNSlices() const; // Philippe was here

private:

    double startField_m;
    std::size_t nSlices_m;

    // Not implemented.
    void operator=(const Undulator &);
};

#endif // CLASSIC_Undulator_HH
