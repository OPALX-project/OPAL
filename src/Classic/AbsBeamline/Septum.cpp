//
// Class Septum
//   Interface for a septum magnet
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
#include "AbsBeamline/Septum.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"

extern Inform *gmsg;

Septum::Septum():Septum("")
{}

Septum::Septum(const std::string &name):
    PluginElement(name),
    width_m(0.0) {
}

Septum::Septum(const Septum &right):
    PluginElement(right),
    width_m(right.width_m) {
    setGeom(width_m);
}

Septum::~Septum() {}

void Septum::accept(BeamlineVisitor &visitor) const {
    visitor.visitSeptum(*this);
}

void Septum::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + 0.005;
    startField -= 0.005;
    initialise(bunch);
}

void Septum::doInitialise(PartBunchBase<double, 3>* /*bunch*/) {
    *gmsg << "Septum initialise" << endl;
}

double Septum::getWidth() const {
    return width_m;
}

void Septum::setWidth(double width) {
    width_m = width;
    setGeom(width_m);
}

bool Septum::doPreCheck(PartBunchBase<double, 3> *bunch) {
    Vector_t rmin;
    Vector_t rmax;
    bunch->get_bounds(rmin, rmax);
    // interested in absolute maximum
    double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
    double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
    double rbunch_max = std::hypot(xmax, ymax);

    if (rbunch_max > rstart_m - 0.1)  {
        return true;
    }
    return false;
}

bool Septum::doCheck(PartBunchBase<double, 3> *bunch, const int /*turnnumber*/, const double /*t*/, const double /*tstep*/) {

    bool flag = false;
    const double slope = (yend_m - ystart_m) / (xend_m - xstart_m);
    const double halfLength = width_m / 2.0 * std::hypot(slope, 1);
    const double intcept = ystart_m - slope * xstart_m;
    const double intcept1 = intcept - halfLength;
    const double intcept2 = intcept + halfLength;

    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
        const Vector_t& R = bunch->R[i];

        double line1 = std::abs(slope * R(0) + intcept1);
        double line2 = std::abs(slope * R(0) + intcept2);

        if (std::abs(R(1)) > line2 &&
           std::abs(R(1)) < line1 &&
           R(0) > xstart_m    &&
           R(0) < xend_m      &&
           R(1) > ystart_m    &&
           R(1) < yend_m) {

            lossDs_m->addParticle(R, bunch->P[i], bunch->ID[i]);
            bunch->Bin[i] = -1;
            flag = true;
        }
    }
    return flag;
}

ElementBase::ElementType Septum::getType() const {
    return SEPTUM;
}