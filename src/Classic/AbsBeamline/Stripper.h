#ifndef CLASSIC_Stripper_HH
#define CLASSIC_Stripper_HH

#include "AbsBeamline/PluginElement.h"

// Class Stripper
// ------------------------------------------------------------------------
//  Class Stripper defines the abstract interface for a striping foil.

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
