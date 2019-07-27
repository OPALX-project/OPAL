#ifndef CLASSIC_Monitor_HH
#define CLASSIC_Monitor_HH

// ------------------------------------------------------------------------
// $RCSfile: Monitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Monitor
//   Defines the abstract interface for a beam position monitor.
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
#include "BeamlineGeometry/StraightGeometry.h"
#include "Structure/LossDataSink.h"

#include <map>
#include <string>

template <class T, unsigned Dim>
class PartBunchBase;
class BeamlineVisitor;

// Class Monitor
// ------------------------------------------------------------------------
/// Interface for beam position monitors.
//  Class Monitor defines the abstract interface for general beam position
//  monitors.

class Monitor: public Component {

public:

    /// Plane selection.
    enum Plane {
        /// Monitor is off (inactive).
        OFF,
        /// Monitor acts on x-plane.
        X,
        /// Monitor acts on y-plane.
        Y,
        /// Monitor acts on both planes.
        XY
    };

    enum Type {
        TEMPORAL,
        SPATIAL
    };

    /// Constructor with given name.
    explicit Monitor(const std::string &name);

    Monitor();
    Monitor(const Monitor &);
    virtual ~Monitor();

    /// Apply visitor to Monitor.
    virtual void accept(BeamlineVisitor &) const override;

    /// Get geometry.
    virtual StraightGeometry &getGeometry() override = 0;

    /// Get geometry. Version for const object.
    virtual const StraightGeometry &getGeometry() const override = 0;

    /// Get plane on which monitor observes.
    virtual Plane getPlane() const = 0;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool applyToReferenceParticle(const Vector_t &R,
                                          const Vector_t &P,
                                          const double &t,
                                          Vector_t &E,
                                          Vector_t &B) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double &kineticEnergy) override;

    virtual void goOffline() override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    void setOutputFN(std::string fn);

    void setType(Type type);

    static void writeStatistics();
private:

    // Not implemented.
    void operator=(const Monitor &);
    std::string filename_m;               /**< The name of the outputfile*/
    Plane plane_m;
    Type type_m;
    unsigned int numPassages_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    static std::map<double, SetStatistics> statFileEntries_sm;
    static const double halfLength_s;
};

inline
void Monitor::setType(Monitor::Type type) {
    type_m = type;
}

#endif // CLASSIC_Monitor_HH