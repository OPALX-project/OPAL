//
// Class Septum
//   Interface for a septum magnet
//
// Copyright (c) 2016-2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_Septum_HH
#define CLASSIC_Septum_HH

#include "AbsBeamline/PluginElement.h"

class Septum: public PluginElement {

public:
    /// Constructor with given name.
    explicit Septum(const std::string &name);

    Septum();
    Septum(const Septum &);
    void operator=(const Septum &) = delete;
    virtual ~Septum();

    /// Apply visitor to Septum.
    virtual void accept(BeamlineVisitor &) const override;
    ///@{ Override implementation of PluginElement
    virtual ElementBase::ElementType getType() const override;
    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;
    ///@}
    /// unhide PluginElement::initialise(PartBunchBase<double, 3> *bunch)
    using PluginElement::initialise;

    ///@{ Member variable access
    void   setWidth(double width);
    double getWidth() const;
    ///@}

private:
    /// Hook for initialise
    virtual void doInitialise(PartBunchBase<double, 3> *bunch) override;
    /// Record hits when bunch particles pass
    virtual bool doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) override;
    /// Virtual hook for preCheck
    virtual bool doPreCheck(PartBunchBase<double, 3>*) override;

    ///@{ input geometry positions
    double width_m;
    ///@}
};

#endif // CLASSIC_Septum_HH
