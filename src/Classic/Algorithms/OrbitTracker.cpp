// ------------------------------------------------------------------------
// $RCSfile: OrbitTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OrbitTracker
//   An visitor class allowing to tracking an orbit through a beamline.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:01:13 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Algorithms/OrbitTracker.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Algorithms/OpalParticle.h"
#include "Algorithms/MapIntegrator.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "Physics/Physics.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Undulator.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Patch.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"

#include <cmath>

typedef FTps<double, 2> Series2;

using Physics::c;


//: Abstract tracker class.
//  An visitor class implementing tracking of an orbit through a beam line.

// Class OrbitTracker
// ------------------------------------------------------------------------


OrbitTracker::OrbitTracker(const Beamline &beamline,
                           const PartData &reference,
                           bool backBeam, bool backTrack):
    AbstractTracker(beamline, reference, backBeam, backTrack)
{}


OrbitTracker::~OrbitTracker()
{}


const FVector<double, 6> &OrbitTracker::getOrbit() const {
    return itsOrbit;
}


void OrbitTracker::setOrbit(const FVector<double, 6> orbit) {
    itsOrbit = orbit;
}


void OrbitTracker::visitBeamBeam(const BeamBeam &) {
    // *** MISSING *** Map for beam-beam.
}


void OrbitTracker::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void OrbitTracker::visitDegrader(const Degrader &deg) {
    applyDrift(flip_s * deg.getElementLength());
}

void OrbitTracker::visitParallelPlate(const ParallelPlate &pplate) {
    //do nothing in orbittracker.
}


void OrbitTracker::visitCyclotronValley(const CyclotronValley &cv) {
    //do nothing in orbittracker.
}

void OrbitTracker::visitComponent(const Component &comp) {
    std::unique_ptr<PartBunchBase<double, 3> > bunch(new PartBunch(&itsReference));
    OpalParticle part(itsOrbit[X], itsOrbit[PX],
                      itsOrbit[Y], itsOrbit[PY],
                      itsOrbit[T], itsOrbit[PT]);
    bunch->push_back(part);
    comp.trackBunch(bunch.get(), itsReference, back_beam, back_track);
    itsOrbit[X]  = part.x();
    itsOrbit[PX] = part.px();
    itsOrbit[Y]  = part.y();
    itsOrbit[PY] = part.py();
    itsOrbit[T]  = part.t();
    itsOrbit[PY] = part.pt();
}


void OrbitTracker::visitCorrector(const Corrector &corr) {
    // Drift through first half of length.
    double length = flip_s * corr.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply kick.
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BDipoleField &field = corr.getField();
    itsOrbit[PX] -= field.getBy() * scale;
    itsOrbit[PY] += field.getBx() * scale;

    // Drift through second half of length.
    if(length) applyDrift(length / 2.0);
}


void OrbitTracker::visitDiagnostic(const Diagnostic &diag) {
    applyDrift(flip_s * diag.getElementLength());
}


void OrbitTracker::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void OrbitTracker::visitUndulator(const Undulator &u) {
  INFOMSG("Undulator not implemented");
  // ada    applyDrift(flip_s * drift.getElementLength());
}

void OrbitTracker::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void OrbitTracker::visitLambertson(const Lambertson &lamb) {
    // Assume the orbit goes through the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void OrbitTracker::visitMarker(const Marker &) {
    // Do nothing.
}


void OrbitTracker::visitMonitor(const Monitor &moni) {
    applyDrift(flip_s * moni.getElementLength());
}


void OrbitTracker::visitMultipole(const Multipole &mult) {
    double length = flip_s * mult.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = mult.getField();

    if(length) {
        // Normal case: Finite-length multipole, field coefficients are B.
        applyMultipoleBody(length, length, field, scale);
    } else {
        // Special case: Thin multipole, field coefficients are integral(B*dl).
        scale *= flip_s;
        applyThinMultipole(field, scale);
    }
}


void OrbitTracker::visitPatch(const Patch &patch) {
    Euclid3D transform = patch.getPatch();
    if(back_path) transform = Inverse(transform);
    applyTransform(transform, 0.0);
}

void OrbitTracker::visitProbe(const Probe &prob) {
    // Do nothing.
}

void OrbitTracker::visitRBend(const RBend &bend) {
    const RBendGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();

    if(length == 0.0) {
        double half_angle = flip_s * geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);

        // Transform from in-plane to mid-plane.
        applyTransform(rotat, 0.0);

        // Apply multipole kick.
        applyThinMultipole(field, scale);

        // Transform from mid-plane to out-plane.
        applyTransform(rotat, 0.0);
    } else {
        double refLength = flip_s * geometry.getArcLength();

        if(back_path) {
            // Transform from global to local.
            applyTransform(Inverse(geometry.getExitPatch()), 0.0);

            // Apply entrance fringe field.
            applyEntranceFringe(bend.getExitFaceRotation(), field, scale);

            // Traverse rbend body.
            applyMultipoleBody(length, refLength, field, scale);

            // Apply exit fringe field.
            applyExitFringe(bend.getEntryFaceRotation(), field, scale);

            // Apply rotation local to global.
            applyTransform(Inverse(geometry.getEntrancePatch()), 0.0);
        } else {
            // Apply rotation global to local.
            applyTransform(geometry.getEntrancePatch(), 0.0);

            // Apply entrance fringe field.
            applyEntranceFringe(bend.getEntryFaceRotation(), field, scale);

            // Traverse rbend body.
            applyMultipoleBody(length, refLength, field, scale);

            // Apply exit fringe field.
            applyExitFringe(bend.getExitFaceRotation(), field, scale);

            // Apply rotation local to global.
            applyTransform(geometry.getExitPatch(), 0.0);
        }
    }
}


void OrbitTracker::visitRFCavity(const RFCavity &as) {
    // Drift through half length.
    double length = flip_s * as.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    double time = itsOrbit[T] / (c * itsReference.getBeta());
    double peak = flip_s * as.getAmplitude() / itsReference.getP();
    double phase = as.getPhase() + as.getFrequency() * time;
    itsOrbit[PT] += peak * sin(phase);

    // Drift through half length.
    if(length) applyDrift(length / 2.0);
}


void OrbitTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Map for RF Quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void OrbitTracker::visitSBend(const SBend &bend) {
    const PlanarArcGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();

    if(length == 0.0) {
        double half_angle = geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);

        // Transform from in-plane to mid-plane.
        applyTransform(rotat, 0.0);

        // Apply multipole kick.
        // Curvature is unknown for zero-length bend.
        applyThinMultipole(field, scale);

        // Transform from mid-plane to out-plane.
        applyTransform(rotat, 0.0);
    } else {
        double h = geometry.getCurvature();
        double refLength = flip_s * geometry.getArcLength();

        if(back_path) {
            // Apply entrance fringe field.
            applyEntranceFringe(bend.getExitFaceRotation(), field, scale);

            // Traverse sbend body.
            applySBendBody(length, refLength, h, field, scale);

            // Apply exit fringe field.
            applyExitFringe(bend.getEntryFaceRotation(), field, scale);
        } else {
            // Apply entrance fringe field.
            applyEntranceFringe(bend.getEntryFaceRotation(), field, scale);

            // Traverse sbend body.
            applySBendBody(length, refLength, h, field, scale);

            // Apply exit fringe field.
            applyExitFringe(bend.getExitFaceRotation(), field, scale);
        }
    }
}


void OrbitTracker::visitSeparator(const Separator &sep) {
    // Drift through first half of length.
    double length = flip_s * sep.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Electrostatic kick.
    double scale = (length * itsReference.getQ()) / itsReference.getP();
    double Ex = scale * sep.getEx();
    double Ey = scale * sep.getEy();
    double pt = 1.0 + itsOrbit[PT];
    itsOrbit[PX] += Ex / pt;
    itsOrbit[PY] += Ey / pt;

    if(length) applyDrift(length / 2.0);
}


void OrbitTracker::visitSeptum(const Septum &sept) {
    // Assume the orbit goes through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void OrbitTracker::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double C = cos(ks * length);
            double S = sin(ks * length);

            double xt  = C * itsOrbit[X]  + S * itsOrbit[Y];
            double yt  = C * itsOrbit[Y]  - S * itsOrbit[X];
            double pxt = C * itsOrbit[PX] + S * itsOrbit[PY];
            double pyt = C * itsOrbit[PY] - S * itsOrbit[PX];

            itsOrbit[X]  = C * xt  + (S / ks) * pxt;
            itsOrbit[Y]  = C * yt  + (S / ks) * pyt;
            itsOrbit[PX] = C * pxt - (S * ks) * xt;
            itsOrbit[PY] = C * pyt - (S * ks) * yt;

            double kin = itsReference.getM() / itsReference.getP();
            itsOrbit[T] += length * itsOrbit[PT] * kin * kin;
        } else {
            applyDrift(length);
        }
    }
}


void OrbitTracker::visitAlignWrapper(const AlignWrapper &wrap) {
    if(wrap.offset().isIdentity()) {
        wrap.getElement()->accept(*this);
    } else {
        Euclid3D e1 = wrap.getEntranceTransform();
        Euclid3D e2 = wrap.getExitTransform();

        if(back_path) {
            // Tracking right to left.
            applyTransform(Inverse(e2), 0.0);
            wrap.getElement()->accept(*this);
            applyTransform(Inverse(e1), 0.0);
        } else {
            // Tracking left to right.
            applyTransform(e1, 0.0);
            wrap.getElement()->accept(*this);
            applyTransform(e2, 0.0);
        }
    }
}


void OrbitTracker::applyDrift(double length) {
    double kin = itsReference.getM() / itsReference.getP();
    itsOrbit[X] += length * itsOrbit[PX];
    itsOrbit[Y] += length * itsOrbit[PY];
    itsOrbit[T] += length * itsOrbit[PT] * kin * kin;
}


void OrbitTracker::applyEntranceFringe
(double edge, const BMultipoleField &field, double scale) {
    double hx = scale * field.normal(1);
    double ex = hx * tan(edge);
    double ey = hx * tan(edge + itsOrbit[PX]);
    itsOrbit[PX] += ex * itsOrbit[X];
    itsOrbit[PY] -= ey * itsOrbit[Y];
}


void OrbitTracker::applyExitFringe
(double edge, const BMultipoleField &field, double scale) {
    double hx = scale * field.normal(1);
    double ex = hx * tan(edge);
    double ey = hx * tan(edge - itsOrbit[PX]);
    itsOrbit[PX] += ex * itsOrbit[X];
    itsOrbit[PY] -= ey * itsOrbit[Y];
}


void OrbitTracker::applyLinearMap
(double length, double refLength, double h,
 const FTps<double, 2> &Fx, const FTps<double, 2> &Fy) {
    // Extract the phase coordinates.
    double x  = itsOrbit[X];
    double px = itsOrbit[PX];
    double y  = itsOrbit[Y];
    double py = itsOrbit[PY];
    double pt = itsOrbit[PT];
    double pt1 = 1.0 - pt;

    // Extract coefficients for equations of motion.
    // These coefficients are w. r. t. the position (x,y).
    // Indexing: 0 = constant, 1 = X, 2 = Y.
    double kx = Fx[1] * pt1;
    double ks = (Fx[2] - Fy[1]) * (pt1 / 2.0);
    double ky = - Fy[2] * pt1;
    double hx = h - Fx[0] * pt1;
    double hy = Fy[0] * pt1;

    // Parameters for longitudinal motion.
    double kin = itsReference.getM() / itsReference.getE();
    double refTime = refLength * kin * kin;

    // Test for zero skew quadrupole component.
    if(ks == 0.0) {
        // Transport coefficients.
        double cx, sx, dx, fx, cy, sy, dy, fy;
        makeFocus(kx, length, cx, sx, dx, fx);
        makeFocus(ky, length, cy, sy, dy, fy);

        // Advance through field.
        itsOrbit[X]  = sx * px + dx * hx;
        itsOrbit[PX] = cx * px + sx * hx;
        itsOrbit[Y]  = sy * py + dy * hy;
        itsOrbit[PY] = cy * py + sy * hy;
        itsOrbit[T] += h * (dx * px + fx * hx);
    } else {
        // Find transformation to principal axes.
        double s1 = (kx + ky) / 2.0;
        double d1 = (kx - ky) / 2.0;
        double root = sqrt(d1 * d1 + ks * ks);
        double c2 = d1 / root;
        double s2 = ks / root;

        // Transport coefficients.
        double cu, su, du, fu, cv, sv, dv, fv;
        double ku = s1 + (d1 * d1 - ks * ks) / root;
        double kv = s1 - (d1 * d1 - ks * ks) / root;
        makeFocus(ku, length, cu, su, du, fu);
        makeFocus(kv, length, cv, sv, dv, fv);

        // Rotate the coordinates to orientation of quadrupole.
        double pu = c2 * px - s2 * py;
        double pv = c2 * py + s2 * px;
        double hu = c2 * (h + hx) - s2 * hy;
        double hv = c2 * hy + s2 * (h + hx);

        // Advance through field.
        itsOrbit[X]  = ((su + sv) * px + (su - sv) * pu +
                        (du + dv) * hx + (du - dv) * hu) / 2.0;
        itsOrbit[PX] = ((cu + cv) * px + (cu - cv) * pu +
                        (su + sv) * hx + (su - sv) * hu) / 2.0;
        itsOrbit[Y]  = ((su + sv) * py - (su - sv) * pv +
                        (du + dv) * hy - (du - dv) * hv) / 2.0;
        itsOrbit[PY] = ((cu + cv) * py - (cu - cv) * pv +
                        (su + sv) * hy - (su - sv) * hv) / 2.0;
        itsOrbit[T] += ((du + dv) * px + (du - dv) * pu +
                        (fu + fv) * hx + (fu - fv) * hu) * (h / 2.0);
    }

    // Add in constant terms.
    itsOrbit[X] += x;
    itsOrbit[Y] += y;
    itsOrbit[T] += refTime * pt + length * h * x;
}


void OrbitTracker::applyMultipoleBody(double length, double refLength,
                                      const BMultipoleField &field, double scale) {
    // Determine normalised field coefficients around actual orbit.
    // Fx and Fy are the normalised field coefficients,
    // expressed as linear functions in x and y.
    int order = field.order();

    if(order > 0) {
        static Series2 x = Series2::makeVariable(0);
        static Series2 y = Series2::makeVariable(1);

        Series2 Fx =   field.normal(order);
        Series2 Fy = - field.skew(order);

        while(order > 1) {
            Series2 Fxt = x * Fx - y * Fy;
            Series2 Fyt = x * Fy + y * Fx;
            --order;
            Fx = Fxt + field.normal(order);
            Fy = Fyt - field.skew(order);
        };
        Fx.setTruncOrder(Fx.getMaxOrder());
        Fy.setTruncOrder(Fy.getMaxOrder());

        Fx *= scale;
        Fy *= scale;

        FVps<double, 2> to_fixpt;
        to_fixpt[0] += itsOrbit[X];
        to_fixpt[1] += itsOrbit[Y];

        Fx = Fx.substitute(to_fixpt);
        Fy = Fy.substitute(to_fixpt);
        applyLinearMap(length, refLength, 0.0, Fx, Fy);
    } else applyDrift(length);
}


void OrbitTracker::applySBendBody(double length, double refLength, double h,
                                  const BMultipoleField &field, double scale) {
    // Determine normalised field coefficients around actual orbit.
    // As is the vector potential times (1 + h*x),
    // expressed as a linear function in x and y.
    Series2 As = buildSBendVectorPotential(field, h) * scale;

    // Fx and Fy are the normalised field coefficients times (1 + h*x).
    Series2 Fx = + As.derivative(0);
    Series2 Fy = - As.derivative(1);
    applyLinearMap(length, refLength, h, Fx, Fy);
}


void OrbitTracker::applyThinMultipole
(const BMultipoleField &field, double factor) {
    int order = field.order();

    if(order > 0) {
        double x = itsOrbit[X];
        double y = itsOrbit[Y];
        double kx = + field.normal(order);
        double ky = - field.skew(order);

        while(--order > 0) {
            double kxt = x * kx - y * ky;
            double kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        itsOrbit[PX] -= kx * factor;
        itsOrbit[PY] += ky * factor;
    }
}


void OrbitTracker::applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        double px1 = itsOrbit[PX];
        double py1 = itsOrbit[PY];
        double pt  = itsOrbit[PT] + 1.0;
        double pz1 = sqrt(pt * pt - px1 * px1 - py1 * py1);

        itsOrbit[PX] = euclid.M(0, 0) * px1 + euclid.M(1, 0) * py1 + euclid.M(2, 0) * pz1;
        itsOrbit[PY] = euclid.M(0, 1) * px1 + euclid.M(1, 1) * py1 + euclid.M(2, 1) * pz1;
        double pz2 = euclid.M(0, 2) * px1 + euclid.M(1, 2) * py1 + euclid.M(2, 2) * pz1;

        double x = itsOrbit[X] - euclid.getX();
        double y = itsOrbit[Y] - euclid.getY();
        double x2 =
            euclid.M(0, 0) * x + euclid.M(1, 0) * y + euclid.M(2, 0) * euclid.getZ();
        double y2 =
            euclid.M(0, 1) * x + euclid.M(1, 1) * y + euclid.M(2, 1) * euclid.getZ();
        double s2 =
            euclid.M(0, 2) * x + euclid.M(1, 2) * y + euclid.M(2, 2) * euclid.getZ();
        double sByPz = s2 / pz2;

        double kin = itsReference.getM() / itsReference.getP();
        double E = sqrt(pt * pt + kin * kin);
        double refTime = refLength / itsReference.getBeta();
        itsOrbit[X] = x2 - sByPz * itsOrbit[PX];
        itsOrbit[Y] = y2 - sByPz * itsOrbit[PY];
        itsOrbit[T] += pt * (- refTime / E  + sByPz);
    }
}


Series2
OrbitTracker::buildSBendVectorPotential(const BMultipoleField &field, double h) {
    // Check sanity.
    if(h == 0.) {
        std::cerr << " <*** ERROR ***> in LinearMapper::buildSBendVectorPotential():\n"
                  << "   attempt to use an infinite radius of curvature." << std::endl;
        throw DomainError("buildSBendVectorPotential(const BMultipoleField &, double)");
    }

    int order = field.order();
    Series2 As;

    if(order > 0) {
        static Series2 x = Series2::makeVariable(0);
        static Series2 y = Series2::makeVariable(1);

        // Terms even in y.
        Series2 Ae = + field.normal(order);
        // Terms odd in y.
        Series2 Ao = - field.skew(order);

        int i = order;
        while(i > 1) {
            i--;
            Ae = Ae * x + field.normal(i);
            Ao = Ao * x - field.skew(i);
        };
        Ae.setTruncOrder(Ae.getMaxOrder());
        Ao.setTruncOrder(Ao.getMaxOrder());

        {
            Series2 hx1 = 1. + h * x;
            Ae = + (Ae * hx1).integral(0);
            Ao = - (Ao * hx1);
        }
        // Add terms up to maximum order.
        As = Ae + y * Ao;

        int k = 2;
        if(k <= order) {
            Series2 yp = y * y / 2.0;
            Series2 radius = 1.0 / h + x;

            while(true) {
                // Terms even in y.
                Ae = Ae.derivative(0);
                Ae = Ae / radius - Ae.derivative(0); // replace this line by following if we decide to treat case h==0.
                // Ae = (h==0. ? 0. : Ae/radius) - Ae.derivative(0);
                As += Ae * yp;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(0);
                Ao = Ao / radius - Ao.derivative(0); // replace this line by following if we decide to treat case h==0.
                // Ao = (h==0. ? 0. : Ao/radius) - Ao.derivative(0);
                As += Ao * yp;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    FVps<double, 2> to_fixpt;
    to_fixpt[0] += itsOrbit[X];
    to_fixpt[1] += itsOrbit[Y];

    return As.substitute(to_fixpt);

}


void OrbitTracker::makeFocus
(double k, double L, double &c, double &s, double &d, double &f) {
    double t = k * L * L;
    if(std::abs(t) < 1.0e-4) {
        c = 1.0 - t / 2.0;
        s = L * (1.0 - t / 6.0);
        d = L * L * (0.5 - t / 24.0);
        f = L * L * L * ((1.0 / 6.0) - t / 120.0);
    } else if(k > 0.0) {
        double r = sqrt(k);
        c = cos(r * L);
        s = sin(r * L) / r;
        d = (1.0 - c) / k;
        f = (L - s) / k;
    } else {
        double r = sqrt(- k);
        c = cosh(r * L);
        s = sinh(r * L) / r;
        d = (1.0 - c) / k;
        f = (L - s) / k;
    }
}
