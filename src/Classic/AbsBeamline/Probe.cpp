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
#include "AbsBeamline/Probe.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Structure/PeakFinder.h"

extern Inform *gmsg;

Probe::Probe():Probe("")
{}

Probe::Probe(const std::string &name):
    PluginElement(name),
    step_m(0.0)
{}

Probe::Probe(const Probe &right):
    PluginElement(right),
    step_m(right.step_m)
{}

Probe::~Probe() {}

void Probe::accept(BeamlineVisitor &visitor) const {
    visitor.visitProbe(*this);
}

void Probe::doInitialise(PartBunchBase<double, 3> *bunch) {
    *gmsg << "* Initialize probe" << endl;
    bool singlemode = (bunch->getTotalNum() == 1) ? true : false;
    peakfinder_m = std::unique_ptr<PeakFinder> (new PeakFinder(getOutputFN(), rmin_m, rend_m, step_m, singlemode));
}

void Probe::doGoOffline() {
    *gmsg << "* Probe goes offline " << getName() << endl;
    if (online_m && peakfinder_m)
        peakfinder_m->save();
    peakfinder_m.reset(nullptr);
}

void Probe::setStep(double step) {
    step_m = step;
}

double Probe::getStep() const {
    return step_m;
}

bool Probe::doPreCheck(PartBunchBase<double, 3> *bunch) {
    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);
    // interested in absolute minimum and maximum
    double xmin = std::min(std::abs(rmin(0)), std::abs(rmax(0)));
    double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
    double ymin = std::min(std::abs(rmin(1)), std::abs(rmax(1)));
    double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
    double rbunch_min = std::hypot(xmin, ymin);
    double rbunch_max = std::hypot(xmax, ymax);

    if (rbunch_max > rmin_m - 0.01 && rbunch_min < rend_m + 0.01 ) {
        return true;
    }
    return false;
}

bool Probe::doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {
    Vector_t probepoint;
    size_t tempnum = bunch->getLocalNum();

    for (unsigned int i = 0; i < tempnum; ++i) {
        double tangle = calculateIncidentAngle(bunch->P[i](0), bunch->P[i](1));
        changeWidth(bunch, i, tstep, tangle);
        int pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
        if (pflag == 0) continue;
        // calculate closest point at probe -> better to use momentum direction
        // probe: y = -A/B * x - C/B or A*X + B*Y + C = 0
        // perpendicular line through R: y = B/A * x + R(1) - B/A * R(0)
        // probepoint(0) = (B_m*B_m*bunch->R[i](0) - A_m*B_m*bunch->R[i](1) - A_m*C_m) / (R_m*R_m);
        // probepoint(1) = (A_m*A_m*bunch->R[i](1) - A_m*B_m*bunch->R[i](0) - B_m*C_m) / (R_m*R_m);
        // probepoint(2) = bunch->R[i](2);
        // calculate time correction for probepoint
        // dist1 > 0, right hand, dt > 0; dist1 < 0, left hand, dt < 0
        double dist1 = (A_m * bunch->R[i](0) + B_m * bunch->R[i](1) + C_m) / R_m; // [m]
        double dist2 = dist1 * std::sqrt(1.0 + 1.0 / tangle / tangle);
        double dt = dist2 / (std::sqrt(1.0 - 1.0 / (1.0 + dot(bunch->P[i], bunch->P[i]))) * Physics::c) * 1.0e9;

        probepoint = bunch->R[i] + dist2 * bunch->P[i] / euclidean_norm(bunch->P[i]);

        // peak finder uses millimetre not metre
        peakfinder_m->addParticle(probepoint * 1e3);

        lossDs_m->addParticle(probepoint, bunch->P[i], bunch->ID[i], t+dt,
                              turnnumber, bunch->bunchNum[i]);
    }

    peakfinder_m->evaluate(turnnumber);

    // we do not lose particles in the probe
    return false;
}

ElementBase::ElementType Probe::getType() const {
    return PROBE;
}