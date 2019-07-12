#ifndef CLASSIC_PluginElement_HH
#define CLASSIC_PluginElement_HH

#include "AbsBeamline/Component.h"
#include <string>
#include <memory>

template <class T, unsigned Dim>
class PartBunchBase;

class LossDataSink;

// Class PluginElement
// ------------------------------------------------------------------------
/// Abstract Interface for (Cyclotron) Plugin Elements (Probe, Stripper, Septum etc.)
/// Implementation via Non-Virtual Interface Template Method

class PluginElement: public Component {

public:
    /// Constructor with given name.
    explicit PluginElement(const std::string &name);

    PluginElement();
    PluginElement(const PluginElement &);
    void operator=(const PluginElement &) = delete;
    virtual ~PluginElement();

    ///@{ Pure virtual implementation of Component
    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override; // not used?
    void initialise(PartBunchBase<double, 3> *bunch); // replacement for virtual initialise
    virtual void finalise()  final; // final since virtual hook doFinalise
    virtual void goOffline() final; // final since virtual hook doGoOffline
    virtual bool bends() const override;
    virtual void getDimensions(double &zBegin, double &zEnd) const override;
    ///@}
    ///@{ Virtual implementation of Component
    virtual bool apply(const size_t &i,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B) override;

    virtual bool applyToReferenceParticle(const Vector_t &R,
                                          const Vector_t &P,
                                          const double &t,
                                          Vector_t &E,
                                          Vector_t &B) override;
    ///@}

    /// Set output filename
    void setOutputFN(std::string fn);
    /// Get output filename
    std::string getOutputFN();

    /// Set dimensions and consistency checks
    void setDimensions(double xstart, double xend, double ystart, double yend);

    ///@{ Member variable access
    double getXStart() const;
    double getXEnd()   const;
    double getYStart() const;
    double getYEnd()   const;
    ///@}
    /// Check if bunch particles are lost
    bool check(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep);
    /// Checks if coordinate is within element
    int  checkPoint(const double & x, const double & y) const;

protected:
    /// Sets geometry geom_m with element width dist
    void setGeom(const double dist);
    /// Change probe width depending on step size and angle of particle
    void changeWidth(PartBunchBase<double, 3> *bunch, const double tstep);
    /// Calculate angle of particle/bunch wrt to element
    double calculateIncidentAngle(double xp, double yp) const;

private:
    /// Check if bunch is close to element
    bool preCheck(PartBunchBase<double, 3> *bunch) {return doPreCheck(bunch);}
    /// Finalise call after check
    bool finaliseCheck(PartBunchBase<double, 3> *bunch, bool flagNeedUpdate) {return doFinaliseCheck(bunch, flagNeedUpdate);}
    /// Pure virtual hook for initialise
    virtual void doInitialise(PartBunchBase<double, 3> *bunch) {}
    /// Pure virtual hook for check
    virtual bool doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) = 0;
    /// Virtual hook for setGeom
    virtual void doSetGeom() {};
    /// Virtual hook for preCheck
    virtual bool doPreCheck(PartBunchBase<double, 3>*) {return true;}
    /// Virtual hook for finaliseCheck
    virtual bool doFinaliseCheck(PartBunchBase<double, 3> *, bool flagNeedUpdate) {return flagNeedUpdate;}
    /// Virtual hook for finalise
    virtual void doFinalise() {};
    /// Virtual hook for goOffline
    virtual void doGoOffline() {};

protected:
    /* Members */
    std::string filename_m; /**< The name of the outputfile*/
    double position_m;
    ///@{ input geometry positions
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double rstart_m;
    double rend_m;
    ///@}
    double rmin_m;    ///< radius closest to the origin
    Point  geom_m[5]; ///< actual geometry positions with adaptive width such that each particle hits element once per turn
    double A_m, B_m, R_m, C_m; ///< Geometric lengths used in calculations

    std::unique_ptr<LossDataSink> lossDs_m;   ///< Pointer to Loss instance
};

#endif // CLASSIC_PluginElement_HH
