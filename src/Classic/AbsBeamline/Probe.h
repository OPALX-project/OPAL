#ifndef CLASSIC_Probe_HH
#define CLASSIC_Probe_HH

#include "AbsBeamline/PluginElement.h"

#include <memory>
#include <string>

class PeakFinder;

// Class Probe
// ------------------------------------------------------------------------
/// Interface for probe.
//  Class Probe defines the abstract interface for a probe.

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
