#ifndef CLASSIC_Drift_HH
#define CLASSIC_Drift_HH

// ------------------------------------------------------------------------
// $RCSfile: Drift.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Drift
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



// Class Drift
// ------------------------------------------------------------------------
/// Interface for drift space.
//  Class Drift defines the abstract interface for a drift space.

class Drift: public Component {

public:

    /// Constructor with given name.
    explicit Drift(const std::string &name);

    Drift();
    Drift(const Drift &right);
    virtual ~Drift();

    /// Apply visitor to Drift.
    virtual void accept(BeamlineVisitor &) const override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    //set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices); // Philippe was here

    //set number of slices for map tracking
    std::size_t getNSlices() const; // Philippe was here

private:

    double startField_m;
    std::size_t nSlices_m;

    // Not implemented.
    void operator=(const Drift &);
};

#endif // CLASSIC_Drift_HH
