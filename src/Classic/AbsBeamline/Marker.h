#ifndef CLASSIC_Marker_HH
#define CLASSIC_Marker_HH

// ------------------------------------------------------------------------
// $RCSfile: Marker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Marker
//   Defines the abstract interface for a marker element.
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


// Class Marker
// ------------------------------------------------------------------------
/// Interface for a marker.
//  Class Marker defines the abstract interface for a marker element.

class Marker: public Component {

public:

    /// Constructor with given name.
    explicit Marker(const std::string &name);

    Marker();
    Marker(const Marker &);
    virtual ~Marker();

    /// Apply visitor to Marker.
    virtual void accept(BeamlineVisitor &) const override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

private:

    // Not implemented.
    void operator=(const Marker &);
};

#endif // CLASSIC_Marker_HH
