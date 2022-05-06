//
// Class Monitor
//   Defines the abstract interface for a beam position monitor.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_Monitor_HH
#define CLASSIC_Monitor_HH

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Structure/LossDataSink.h"

#include <map>
#include <string>

template <class T, unsigned Dim>
class PartBunchBase;
class BeamlineVisitor;

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

    virtual ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    void setCollectionType(CollectionType type);

    static void writeStatistics();

    virtual int getRequiredNumberOfTimeSteps() const override;

    virtual bool isInside(const Vector_t &r) const override;
private:

    void driftToCorrectPositionAndSave(const Vector_t& R, const Vector_t& P);

    // Not implemented.
    void operator=(const Monitor &);
    std::string filename_m; /**< The name of the outputfile*/
    Plane plane_m;
    CollectionType type_m;
    unsigned int numPassages_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    static std::map<double, SetStatistics> statFileEntries_sm;
    static const double halfLength_s;
};

inline
void Monitor::setCollectionType(CollectionType type) {
    type_m = type;
}

inline
int Monitor::getRequiredNumberOfTimeSteps() const
{
    return 1;
}

inline
bool Monitor::isInside(const Vector_t &r) const
{
    const double length = getElementLength();
    return std::abs(r(2)) <= 0.5 * length && isInsideTransverse(r);
}


#endif // CLASSIC_Monitor_HH