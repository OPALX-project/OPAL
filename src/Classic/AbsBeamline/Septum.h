#ifndef CLASSIC_Septum_HH
#define CLASSIC_Septum_HH

#include "AbsBeamline/PluginElement.h"

// Class Septum
// ------------------------------------------------------------------------
/// Interface for septum magnet.
//  Class Septum defines the abstract interface for a septum magnet.

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

    ///@{ input geometry positions
    double width_m;
    ///@}
};

#endif // CLASSIC_Septum_HH
