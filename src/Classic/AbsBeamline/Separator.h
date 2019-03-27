#ifndef CLASSIC_Separator_HH
#define CLASSIC_Separator_HH

// ------------------------------------------------------------------------
// $RCSfile: Separator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Separator
//   Defines the abstract interface for an  separator.
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


// Class Separator
// ------------------------------------------------------------------------
/// Interface for electrostatic separator.
//  Class Separator defines the abstract interface for electrostatic
//  separators.

class Separator: public Component {

public:

    /// Constructor with given name.
    explicit Separator(const std::string &name);

    Separator();
    Separator(const Separator &);
    virtual ~Separator();

    /// Apply visitor to Separator.
    virtual void accept(BeamlineVisitor &) const override;

    /// Get horizontal component Ex of field in V/m.
    virtual double getEx() const = 0;

    /// Get vertical component Ey of field in V/m.
    virtual double getEy() const = 0;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

private:

    // Not implemented.
    void operator=(const Separator &);
};

#endif // CLASSIC_Separator_HH
