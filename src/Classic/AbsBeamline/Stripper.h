//
// Class Stripper
//   The Stripper element defines the interface for a stripping foil
//
// Copyright (c) 2011, Jianjun Yang, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2017-2019, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#ifndef CLASSIC_Stripper_HH
#define CLASSIC_Stripper_HH

#include "AbsBeamline/PluginElement.h"

class Stripper: public PluginElement {

public:
    /// Constructor with given name.
    explicit Stripper(const std::string &name);

    Stripper();
    Stripper(const Stripper &);
    void operator=(const Stripper &) = delete;
    virtual ~Stripper();

    /// Apply visitor to Stripper.
    virtual void accept(BeamlineVisitor &) const override;
    ///@{ Override implementation of PluginElement
    virtual ElementBase::ElementType getType() const override;
    ///@}

    ///@{ Member variable access
    void   setOPCharge(double charge);
    double getOPCharge() const;
    void   setOPMass(double mass);
    double getOPMass() const;
    void   setOPYield(double yield);
    double getOPYield() const;
    void   setStop(bool stopflag);
    bool   getStop() const;
    ///@}

private:
    /// Record hits when bunch particles pass
    virtual bool doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) override;
    /// Virtual hook for finalise
    virtual void doFinalise() override;
    /// Virtual hook for preCheck
    virtual bool doPreCheck(PartBunchBase<double, 3>*) override;
    /// Virtual hook for finaliseCheck
    virtual bool doFinaliseCheck(PartBunchBase<double, 3> *bunch, bool flagNeedUpdate) override;

    double opcharge_m; ///< Charge number of the out-coming particle
    double opmass_m;   ///< Mass of the out-coming particle
    double opyield_m;  ///< Yield of the out-coming particle
    bool   stop_m;     ///< Flag if particles should be stripped or stopped
};

#endif // CLASSIC_Stripper_HH
