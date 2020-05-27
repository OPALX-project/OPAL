//
// Class Probe
//   Interface for a probe
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
#ifndef CLASSIC_Probe_HH
#define CLASSIC_Probe_HH

#include "AbsBeamline/PluginElement.h"

#include <memory>
#include <string>

class PeakFinder;

class Probe: public PluginElement {

public:
    /// Constructor with given name.
    explicit Probe(const std::string &name);

    Probe();
    Probe(const Probe &);
    void operator=(const Probe &) = delete;
    virtual ~Probe();

    /// Apply visitor to Probe.
    virtual void accept(BeamlineVisitor &) const override;

    /// Set probe histogram bin width
    void setStep(double step);
    ///@{ Member variable access
    virtual double getStep() const;
    ///@}
    virtual ElementBase::ElementType getType() const override;

private:
    /// Initialise peakfinder file
    virtual void doInitialise(PartBunchBase<double, 3> *bunch) override;
    /// Record probe hits when bunch particles pass
    virtual bool doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) override;
    /// Hook for goOffline
    virtual void doGoOffline() override;
    /// Virtual hook for preCheck
    virtual bool doPreCheck(PartBunchBase<double, 3>*) override;

    double step_m; ///< Step size of the probe (bin width in histogram file)
    std::unique_ptr<PeakFinder> peakfinder_m; ///< Pointer to Peakfinder instance
};

#endif // CLASSIC_Probe_HH
