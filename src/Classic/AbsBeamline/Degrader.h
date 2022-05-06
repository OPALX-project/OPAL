//
// Class Degrader
//   Defines the abstract interface for a beam degrader.
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
#ifndef CLASSIC_Degrader_HH
#define CLASSIC_Degrader_HH

#include "AbsBeamline/Component.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/StraightGeometry.h"

#include <string>
#include <vector>

class Degrader: public Component {

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
    explicit Degrader(const std::string& name);

    Degrader();
    Degrader(const Degrader& rhs);
    virtual ~Degrader();

    /// Apply visitor to Degrader.
    virtual void accept(BeamlineVisitor&) const;

    virtual bool apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B);

    virtual bool applyToReferenceParticle(const Vector_t& R,
                                          const Vector_t& P,
                                          const double& t,
                                          Vector_t& E,
                                          Vector_t& B);

    virtual void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField);

    virtual void initialise(PartBunchBase<double, 3>* bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double& kineticEnergy);

    virtual void goOffline();

    virtual ElementType getType() const;

    virtual void getDimensions(double& zBegin, double& zEnd) const;

    virtual bool isInMaterial(double z);

private:

    // Not implemented.
    void operator=(const Degrader&);

    std::vector<double> PosX_m;
    std::vector<double> PosY_m;
    std::vector<double> PosZ_m;
    std::vector<double> MomentumX_m;
    std::vector<double> MomentumY_m;
    std::vector<double> MomentumZ_m;
    std::vector<double> time_m;
    std::vector<int> id_m;
};

#endif // CLASSIC_Degrader_HH
