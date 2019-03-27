#ifndef CLASSIC_RFQuadrupole_HH
#define CLASSIC_RFQuadrupole_HH

// ------------------------------------------------------------------------
// $RCSfile: RFQuadrupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Class: RFQuadrupole
//   *** MISSING *** RFQuadrupole interface is incomplete.
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class RFQuadrupole
// ------------------------------------------------------------------------
/// Interface for RF Quadrupole.
//  Class RFQuadrupole defines the abstract interface for a RF Quadrupole.

class RFQuadrupole: public Component {

public:

    /// Constructor with given name.
    explicit RFQuadrupole(const std::string &name);

    RFQuadrupole();
    RFQuadrupole(const RFQuadrupole &);
    virtual ~RFQuadrupole();

    /// Apply visitor to RFQuadrupole.
    virtual void accept(BeamlineVisitor &) const override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

private:

    // Not implemented.
    void operator=(const RFQuadrupole &);
};

#endif // CLASSIC_RFQuadrupole_HH
