//
// Class Tracker
//  Track particles or bunches.
//  An abstract base class for all visitors capable of tracking particles
//  through a beam element.
//  [P]
//  Phase space coordinates (in this order):
//  [DL]
//  [DT]x:[DD]
//    horizontal displacement (metres).
//  [DT]p_x/p_r:[DD]
//     horizontal canonical momentum (no dimension).
//  [DT]y:[DD]
//    vertical displacement (metres).
//  [DT]p_y/p_r:[DD]
//    vertical canonical momentum (no dimension).
//  [DT]delta_p/p_r:[DD]
//    relative momentum error (no dimension).
//  [DT]v*delta_t:[DD]
//    time difference delta_t w.r.t. the reference frame which moves with
//    uniform velocity
//  [P]
//    v_r = c*beta_r = p_r/m
//  [P]
//    along the design orbit, multiplied by the instantaneous velocity v of
//    the particle (metres).
//  [/DL]
//  Where
//  [DL]
//  [DT]p_r:[DD]
//    is the constant reference momentum defining the reference frame velocity.
//  [DT]m:[DD]
//    is the rest mass of the particles.
//  [/DL]
//  Other units used:
//  [DL]
//  [DT]reference momentum:[DD]
//    electron-volts.
//  [DT]accelerating voltage:[DD]
//    volts.
//  [DT]separator voltage:[DD]
//    volts.
//  [DT]frequencies:[DD]
//    hertz.
//  [DT]phase lags:[DD]
//    multiples of (2*pi).
//  [/DL]
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Algorithms/Tracker.h"
#include "Fields/BMultipoleField.h"

//FIXME Remove headers and dynamic_cast in readOneBunchFromFile
#include "Algorithms/PartBunch.h"
#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

#include <cfloat>
#include <cmath>
#include <limits>

typedef FTps<double, 2> Series2;
typedef FTps<double, 6> Series;

// Class Tracker
// ------------------------------------------------------------------------


Tracker::Tracker(const Beamline &beamline, const PartData &reference,
                 bool backBeam, bool backTrack):
    Tracker(beamline, nullptr, reference, backBeam, backTrack)
{}


Tracker::Tracker(const Beamline &beamline,
                 PartBunchBase<double, 3> *bunch,
                 const PartData &reference,
                 bool backBeam, bool backTrack):
    AbstractTracker(beamline, reference, backBeam, backTrack),
    itsBeamline_m(beamline),
    itsBunch_m(bunch)
{}


Tracker::~Tracker()
{}


const PartBunchBase<double, 3> *Tracker::getBunch() const {
    return itsBunch_m;
}


void Tracker::addToBunch(const OpalParticle &part) {
    itsBunch_m->push_back(part);
}


//~ void Tracker::setBunch(const PartBunch &bunch) {
    //~ itsBunch_m = &bunch;
//~ }


void Tracker::visitComponent(const Component &comp) {
    comp.trackBunch(itsBunch_m, itsReference, back_beam, back_track);
}


void Tracker::applyDrift(double length) {
    double kin = itsReference.getM() / itsReference.getP();
    double refTime = length / itsReference.getBeta();

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->getParticle(i);
        if(part.getX() != std::numeric_limits<double>::max()) {
            const Vector_t& P = part.getP();
            double lByPz = length / std::sqrt(2 * P[2] * P[2] - dot(P, P));
            part.setX(part.getX() + P[0] * lByPz);
            part.setY(part.getY() + P[1] * lByPz);
            part.setZ(part.getZ() + P[2] * (refTime / std::sqrt(P[2] * P[2] + kin * kin) - lByPz));
        }
        itsBunch_m->setParticle(part, i);
    }
}


void Tracker::applyThinMultipole
(const BMultipoleField &field, double scale) {
    int order = field.order();

    if(order > 0) {
        for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
            OpalParticle part = itsBunch_m->getParticle(i);
            if(part.getX() != std::numeric_limits<double>::max()) {
                double x = part.getX();
                double y = part.getY();
                double kx = + field.normal(order);
                double ky = - field.skew(order);

                int ord = order;
                while(--ord > 0) {
                    double kxt = x * kx - y * ky;
                    double kyt = x * ky + y * kx;
                    kx = kxt + field.normal(ord);
                    ky = kyt - field.skew(ord);
                }
                part.setPx(part.getPx() - kx * scale);
                part.setPy(part.getPy() + ky * scale);
            }
            itsBunch_m->setParticle(part, i);
        }
    }
}


void Tracker::applyThinSBend
(const BMultipoleField &field, double scale, double h) {
    Series2 As = buildSBendVectorPotential2D(field, h) * scale;
    Series2 Fx = As.derivative(0);
    Series2 Fy = As.derivative(1);

    // These substitutions work because As depends on x and y only,
    // and not on px or py.
    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->getParticle(i);
        FVector<double, 2> z;
        z[0] = part.getX();
        z[1] = part.getY();
        part.setPx(part.getPx() - Fx.evaluate(z));
        part.setPy(part.getPy() - Fy.evaluate(z));
        itsBunch_m->setParticle(part, i);
    }
}


void Tracker::applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        double kin = itsReference.getM() / itsReference.getP();
        double refTime = refLength / itsReference.getBeta();

        for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
            OpalParticle part = itsBunch_m->getParticle(i);
            double px = part.getPx();
            double py = part.getPy();
            double pt = part.getPz() + 1.0;
            double pz = std::sqrt(pt * pt - px * px - py * py);

            part.setPx(euclid.M(0, 0) * px + euclid.M(1, 0) * py + euclid.M(2, 0) * pz);
            part.setPy(euclid.M(0, 1) * px + euclid.M(1, 1) * py + euclid.M(2, 1) * pz);
            pz = euclid.M(0, 2) * px + euclid.M(1, 2) * py + euclid.M(2, 2) * pz;

            double x = part.getX() - euclid.getX();
            double y = part.getY() - euclid.getY();
            double x2 =
                euclid.M(0, 0) * x + euclid.M(1, 0) * y - euclid.M(2, 0) * euclid.getZ();
            double y2 =
                euclid.M(0, 1) * x + euclid.M(1, 1) * y - euclid.M(2, 1) * euclid.getZ();
            double s2 =
                euclid.M(0, 2) * x + euclid.M(1, 2) * y - euclid.M(2, 2) * euclid.getZ();
            double sByPz = s2 / pz;

            double E = std::sqrt(pt * pt + kin * kin);
            part.setX(x2 - sByPz * part.getPx());
            part.setY(y2 - sByPz * part.getPy());
            part.setZ(part.getZ() + pt * (refTime / E + sByPz));
            itsBunch_m->setParticle(part, i);
        }
    }
}


Series2 Tracker::
buildMultipoleVectorPotential2D(const BMultipoleField &field) {
    int order = field.order();

    if(order > 0) {
        static const Series2 x = Series2::makeVariable(0);
        static const Series2 y = Series2::makeVariable(1);
        Series2 kx = + field.normal(order) / double(order);
        Series2 ky = - field.skew(order)   / double(order);

        while(order > 1) {
            Series2 kxt = x * kx - y * ky;
            Series2 kyt = x * ky + y * kx;
            order--;
            kx = kxt + field.normal(order) / double(order);
            ky = kyt - field.skew(order)   / double(order);
        }

        Series2 As = x * kx - y * ky;
        As.setTruncOrder(As.getMaxOrder());
        return As;
    } else {
        return Series2(0.0);
    }
}


Series Tracker::
buildMultipoleVectorPotential(const BMultipoleField &field) {
    int order = field.order();

    if(order > 0) {
        static const Series x = Series::makeVariable(X);
        static const Series y = Series::makeVariable(Y);
        Series kx = + field.normal(order) / double(order);
        Series ky = - field.skew(order)   / double(order);

        while(order > 1) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            order--;
            kx = kxt + field.normal(order) / double(order);
            ky = kyt - field.skew(order)   / double(order);
        }

        Series As = x * kx - y * ky;
        As.setTruncOrder(As.getMaxOrder());
        return As;
    } else {
        return Series(0.0);
    }
}


Series2
Tracker::buildSBendVectorPotential2D(const BMultipoleField &field, double h) {
    int order = field.order();
    Series2 As;

    if(order > 0) {
        static const Series2 x = Series2::makeVariable(0);
        static const Series2 y = Series2::makeVariable(1);

        // Construct terms constant and linear in y.
        Series2 Ae = + field.normal(order); // Term even in y.
        Series2 Ao = - field.skew(order);   // Term odd  in y.

        for(int i = order; --i >= 1;) {
            Ae = Ae * x + field.normal(i);
            Ao = Ao * x - field.skew(i);
        }
        Ae.setTruncOrder(Ae.getMaxOrder());
        Ao.setTruncOrder(Ao.getMaxOrder());

        Series2 hx1 = 1. + h * x; // normalized radius
        Ae = + (Ae * hx1).integral(X);
        Ao = - (Ao * hx1);
        // Add terms up to maximum order.
        As = Ae + y * Ao;

        int k = 2;
        if(k <= order) {
            Series2 yp = y * y / 2.0;

            while(true) {
                // Terms even in y.
                Ae = Ae.derivative(0);
                Ae = h * Ae / hx1 - Ae.derivative(0);
                As += Ae * yp;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(0);
                Ao = h * Ao / hx1 - Ao.derivative(0);
                As += Ao * yp;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    return As;
}


Series
Tracker::buildSBendVectorPotential(const BMultipoleField &field, double h) {
    int order = field.order();
    Series As;

    if(order > 0) {
        static const Series x = Series::makeVariable(X);
        static const Series y = Series::makeVariable(Y);

        // Construct terms constant and linear in y.
        Series Ae = + field.normal(order); // Term even in y.
        Series Ao = - field.skew(order);   // Term odd  in y.

        for(int i = order; --i >= 1;) {
            Ae = Ae * x + field.normal(i);
            Ao = Ao * x - field.skew(i);
        }
        Ae.setTruncOrder(Ae.getMaxOrder());
        Ao.setTruncOrder(Ao.getMaxOrder());

        Series hx1 = 1. + h * x; // normalized radius
        Ae = + (Ae * hx1).integral(X);
        Ao = - (Ao * hx1);
        // Add terms up to maximum order.
        As = Ae + y * Ao;

        int k = 2;
        if(k <= order) {
            Series yp = y * y / 2.0;

            while(true) {
                // Terms even in y.
                Ae = Ae.derivative(X);
                Ae = h * Ae / hx1 - Ae.derivative(X);
                As += Ae * yp;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(X);
                Ao = h * Ao / hx1 - Ao.derivative(X);
                As += Ao * yp;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    return As;
}