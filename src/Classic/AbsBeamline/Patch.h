#ifndef CLASSIC_Patch_HH
#define CLASSIC_Patch_HH

// ------------------------------------------------------------------------
// $RCSfile: Patch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Patch
//   Defines the abstract interface for a geometry patch.
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


// Class Patch
// ------------------------------------------------------------------------
/// Interface for a geometric patch.
//  Class Patch defines the abstract interface for geometric patches,
//  i.e. a geometric transform representing a misalignment of local
//  coordinate systems of two subsequent element with each other.

class Patch: public Component {

public:

    /// Constructor with given name.
    explicit Patch(const std::string &name);

    Patch();
    Patch(const Patch &);
    virtual ~Patch();

    /// Apply visitor to patch.
    virtual void accept(BeamlineVisitor &) const override;

    /// Get patch transform.
    virtual const Euclid3D &getPatch() const = 0;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual  ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

private:

    // Not implemented.
    void operator=(const Patch &);
};

#endif // CLASSIC_Patch_HH
