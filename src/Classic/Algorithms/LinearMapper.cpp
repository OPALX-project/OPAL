// ------------------------------------------------------------------------
// $RCSfile: LinearMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.5.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LinearMapper
//   The visitor class for building a linear map for a beamline using
//   linear maps for all elements.
//
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:01:13 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Algorithms/LinearMapper.h"

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

#include "Algorithms/MapIntegrator.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/LinearFun.h"
#include "FixedAlgebra/LinearMath.h"
#include "Physics/Physics.h"

using Physics::c;

typedef FTps<double, 2> Series2;
typedef LinearFun<double, 6> Linear;


// Class LinearMapper
// ------------------------------------------------------------------------

LinearMapper::LinearMapper(const Beamline &beamline, const PartData &reference,
                           bool revBeam, bool revTrack):
    AbstractMapper(beamline, reference, revBeam, revTrack)
{}


LinearMapper::~LinearMapper()
{}


void LinearMapper::getMap(LinearMap<double, 6> &map) const {
    map = itsMap;
}


void LinearMapper::getMap(FVps<double, 6> &map) const {
    map = FVps<double, 6>(itsMap);
}


void LinearMapper::setMap(const LinearMap<double, 6> &map) {
    itsMap = map;
}


void LinearMapper::setMap(const FVps<double, 6> &map) {
    itsMap = LinearMap<double, 6>(map);
}


void LinearMapper::visitBeamBeam(const BeamBeam &) {
    // *** MISSING *** Map for beam-beam.
}

void LinearMapper::visitBeamStripping(const BeamStripping &) {
    // *** MISSING *** Map for beam stripping.
}

void LinearMapper::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}


void LinearMapper::visitComponent(const Component &comp) {
    FVps<double, 6> map(itsMap);
    comp.trackMap(map, itsReference, back_beam, back_track);
    itsMap = LinearMap<double, 6>(map);
}


void LinearMapper::visitCorrector(const Corrector &corr) {
    // Drift through first half of length.
    double length = flip_s * corr.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply kick.
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BDipoleField &field = corr.getField();
    itsMap[PX] -= field.getBy() * scale;
    itsMap[PY] += field.getBx() * scale;

    // Drift through second half of length.
    if(length) applyDrift(length / 2.0);
}

void LinearMapper::visitDegrader(const Degrader &deg) {
    applyDrift(flip_s * deg.getElementLength());
}

void LinearMapper::visitDiagnostic(const Diagnostic &diag) {
    applyDrift(flip_s * diag.getElementLength());
}


void LinearMapper::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void LinearMapper::visitUndulator(/* const Undulator &u */) {
    INFOMSG("Undulator not implemented");
}

void LinearMapper::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void LinearMapper::visitLambertson(const Lambertson &lamb) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void LinearMapper::visitMarker(const Marker &/*marker*/) {
    // Do nothing.
}


void LinearMapper::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void LinearMapper::visitMultipole(const Multipole &multipole) {
    double length = flip_s * multipole.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = multipole.getField();

    if(length) {
        // Normal case: Finite-length multipole, field coefficients are B.
        applyMultipoleBody(length, length, field, scale);
    } else {
        // Special case: Thin multipole, field coefficients are integral(B*dl).
        scale *= flip_s;
        applyThinMultipole(field, scale);
    }
}


void LinearMapper::visitPatch(const Patch &patch) {
    Euclid3D transform = patch.getPatch();
    if(back_track) transform = Inverse(transform);
    applyTransform(transform, 0.0);
}

void LinearMapper::visitProbe(const Probe &/*Prob*/) {
    // Do nothing.
}


#include "Algorithms/rbendmap.h"

void LinearMapper::visitRBend(const RBend &bend) {
    visitRBend0(bend);
}

void LinearMapper::visitRBend0(const RBend &bend) {

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
        scale *= flip_B;
        if(back_track) {
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


void LinearMapper::visitRFCavity(const RFCavity &as) {
    // Drift through half length.
    double length = flip_s * as.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    double freq = as.getFrequency();
    double peak = flip_s * as.getAmplitude() / itsReference.getP();
    Linear phase = as.getPhase() + (freq / c) * itsMap[T];
    itsMap[PT] += peak * sin(phase);

    // Drift through half length.
    if(length) applyDrift(length / 2.0);
}


void LinearMapper::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Map for RF Quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void LinearMapper::visitSBend(const SBend &bend) {
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

        if(back_track) {
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


void LinearMapper::visitSeparator(const Separator &sep) {
    // Drift through first half of length.
    double length = flip_s * sep.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Electrostatic kick.
    double scale = (length * itsReference.getQ()) / itsReference.getP();
    double Ex = scale * sep.getEx();
    double Ey = scale * sep.getEy();
    Linear pt = 1.0 + itsMap[PT];
    itsMap[PX] += Ex / pt;
    itsMap[PY] += Ey / pt;

    if(length) applyDrift(length / 2.0);
}


void LinearMapper::visitSeptum(const Septum &sept) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void LinearMapper::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double C = cos(ks * length);
            double S = sin(ks * length);

            Linear xt  = C * itsMap[X]  + S * itsMap[Y];
            Linear yt  = C * itsMap[Y]  - S * itsMap[X];
            Linear pxt = C * itsMap[PX] + S * itsMap[PY];
            Linear pyt = C * itsMap[PY] - S * itsMap[PX];

            itsMap[X]  = C * xt  + (S / ks) * pxt;
            itsMap[Y]  = C * yt  + (S / ks) * pyt;
            itsMap[PX] = C * pxt - (S * ks) * xt;
            itsMap[PY] = C * pyt - (S * ks) * yt;

            double kin = itsReference.getM() / itsReference.getP();
            itsMap[T] += length * itsMap[PT] * kin * kin;
        } else {
            applyDrift(length);
        }
    }
}


void LinearMapper::visitParallelPlate(const ParallelPlate &/*pplate*/) {
    //do nothing
}


void LinearMapper::visitAlignWrapper(const AlignWrapper &wrap) {
    if(wrap.offset().isIdentity()) {
        wrap.getElement()->accept(*this);
    } else {
        Euclid3D e1 = wrap.getEntranceTransform();
        Euclid3D e2 = wrap.getExitTransform();

        if(back_track) {
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


void LinearMapper::visitMapIntegrator(const MapIntegrator &i) {
    FVps<double, 6> map(itsMap);
    i.trackMap(map, itsReference, back_beam, back_track);
    itsMap = LinearMap<double, 6>(map);
}


void LinearMapper::makeFocus
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


void LinearMapper::applyDrift(double length) {
    double kin = itsReference.getM() / itsReference.getP();
    itsMap[X] += length * itsMap[PX];
    itsMap[Y] += length * itsMap[PY];
    itsMap[T] += length * itsMap[PT] * kin * kin;
}


void LinearMapper::applyEntranceFringe(double angle,
                                       const BMultipoleField &field,
                                       double scale) {
    double hx = scale * field.normal(1);
    double ex = hx * tan(angle);
    double ey = hx * tan(angle + itsMap[PX][0]);
    itsMap[PX] += ex * itsMap[X];
    itsMap[PY] -= ey * itsMap[Y];
}


void LinearMapper::applyExitFringe(double angle,
                                   const BMultipoleField &field,
                                   double scale) {
    double hx = scale * field.normal(1);
    double ex = hx * tan(angle);
    double ey = hx * tan(angle - itsMap[PX][0]);
    itsMap[PX] += ex * itsMap[X];
    itsMap[PY] -= ey * itsMap[Y];
}


void LinearMapper::applyLinearMap(double length, double refLength, double h,
                                  const Series2 &Fx, const Series2 &Fy) {
    // Extract the phase coordinates.
    Linear x  = itsMap[X];
    Linear px = itsMap[PX];
    Linear y  = itsMap[Y];
    Linear py = itsMap[PY];
    Linear pt = itsMap[PT];

    // Split position variables into constant and variable terms.
    double x0 = x[0];
    double y0 = y[0];
    x[0] = 0.0;
    y[0] = 0.0;

    // Extract coefficients for equations of motion.
    // Indexing: 0 = constant, 1 = X, 2 = Y.
    double kx = Fx[1];
    double ks = (Fx[2] - Fy[1]) / 2.0;
    double ky = - Fy[2];
    Linear hx = h * (1.0 + pt) - Fx[0];
    double hy = Fy[0];

    // Parameters for longitudinal motion.
    double kin = itsReference.getM() / itsReference.getE();
    double refTime = refLength * kin * kin;

    // Test for zero skew quadrupole component.
    if(ks == 0.0) {
        // Transport coefficients.
        double cx, sx, dx, fx, cy, sy, dy, fy;
        makeFocus(kx, length, cx, sx, dx, fx);
        makeFocus(ky, length, cy, sy, dy, fy);
        double wx = - kx * sx;
        double wy = - ky * sy;

        // Advance through field.
        itsMap[X]  = cx * x + sx * px + dx * hx;
        itsMap[PX] = wx * x + cx * px + sx * hx;
        itsMap[Y]  = cy * y + sy * py + dy * hy;
        itsMap[PY] = wy * y + cy * py + sy * hy;
        itsMap[T] += h * (sx * x + dx * px + fx * hx);
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
        double wu = - ku * su;
        double wv = - kv * sv;

        // Rotate the coordinates to orientation of quadrupole.
        Linear u  = c2 * x  - s2 * y;
        Linear v  = c2 * y  + s2 * x;
        Linear pu = c2 * px - s2 * py;
        Linear pv = c2 * py + s2 * px;
        Linear hu = c2 * (h + hx) - s2 * hy;
        Linear hv = c2 * hy + s2 * (h + hx);

        // Advance through field.
        itsMap[X]  = ((cu + cv) * x  + (cu - cv) * u  +
                      (su + sv) * px + (su - sv) * pu +
                      (du + dv) * hx + (du - dv) * hu) / 2.0;
        itsMap[PX] = ((wu + wv) * x  + (wu - wv) * u  +
                      (cu + cv) * px + (cu - cv) * pu +
                      (su + sv) * hx + (su - sv) * hu) / 2.0;
        itsMap[Y]  = ((cu + cv) * y  - (cu - cv) * v  +
                      (su + sv) * py - (su - sv) * pv +
                      (du + dv) * hy - (du - dv) * hv) / 2.0;
        itsMap[PY] = ((wu + wv) * y  - (wu - wv) * v  +
                      (cu + cv) * py - (cu - cv) * pv +
                      (su + sv) * hy - (su - sv) * hv) / 2.0;
        itsMap[T] += ((su + sv) * x  + (su - sv) * u  +
                      (du + dv) * px + (du - dv) * pu +
                      (fu + fv) * hx + (fu - fv) * hu) * (h / 2.0);
    }

    // Add in constant terms.
    itsMap[X] += x0;
    itsMap[Y] += y0;
    itsMap[T] += refTime * pt + length * h * x0;
}


Series2
LinearMapper::buildSBendVectorPotential(const BMultipoleField &field, double h) {
    //std::cerr << "==> In buildSBendVectorPotential(const BMultipoleField &field, double h)..."
    //          << std::endl;
    // Check sanity.
    if(h == 0.) {
        std::cerr << " <*** ERROR ***> in LinearMapper::buildSBendVectorPotential():\n"
                  << "   attempt to use an infinite radius of curvature." << std::endl;
        throw DomainError("buildSBendVectorPotential(const BMultipoleField &, double)");
    }

    int order = field.order();
    Series2 As;

    //std::cerr << "order = " << order << std::endl;
    //for (int m = 1; m <= order; ++m) {
    //std::cerr << "Order " << m << ": "
    //          << field.normal(m) << " " << field.skew(m) << std::endl;
    //}

    if(order > 0) {
        static Series2 x = Series2::makeVariable(0);
        static Series2 y = Series2::makeVariable(1);

        // Terms even/odd in y.
        Series2 Ae = + field.normal(order);
        Series2 Ao = - field.skew(order);

        int i = order;
        while(i > 1) {
            --i;
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
                Ae = Ae / radius - Ae.derivative(0); // use following line to include case h==0.
                // Ae = (h==0. ? 0. : Ae/radius) - Ae.derivative(0);
                As += Ae * yp;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(0);
                Ao = Ao / radius - Ao.derivative(0); // use following line to include case h==0.
                // Ao = (h==0. ? 0. : Ao/radius) - Ao.derivative(0);
                As += Ao * yp;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    //std::cerr << " As = " << As << std::endl;
    double x0 = itsMap[X][0], y0 = itsMap[Y][0];
    if(x0 != 0. || y0 != 0.) {
        FVps<double, 2> to_fixpt;
        to_fixpt[0] += x0;
        to_fixpt[1] += y0;
        As = As.substitute(to_fixpt);
    }
    //std::cerr << " As(r-fp) = " << As << std::endl;

    return As;
}


void LinearMapper::applyMultipoleBody(double length, double refLength,
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
        to_fixpt[0] += itsMap[X][0];
        to_fixpt[1] += itsMap[Y][0];

        Fx = Fx.substitute(to_fixpt);
        Fy = Fy.substitute(to_fixpt);
        applyLinearMap(length, refLength, 0.0, Fx, Fy);
    } else applyDrift(length);
}


void LinearMapper::applySBendBody(double length, double refLength, double h,
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


void LinearMapper::applyThinMultipole
(const BMultipoleField &field, double scale) {
    int order = field.order();

    if(order > 0) {
        Linear x = itsMap[X];
        Linear y = itsMap[Y];
        Linear kx = + field.normal(order);
        Linear ky = - field.skew(order);

        while(--order > 0) {
            Linear kxt = x * kx - y * ky;
            Linear kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        itsMap[PX] -= kx * scale;
        itsMap[PY] += ky * scale;
    }
}


void LinearMapper::applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        Linear px1 = itsMap[PX];
        Linear py1 = itsMap[PY];
        Linear pt  = itsMap[PT] + 1.0;
        Linear pz1 = sqrt(pt * pt - px1 * px1 - py1 * py1);

        itsMap[PX] = euclid.M(0, 0) * px1 + euclid.M(1, 0) * py1 + euclid.M(2, 0) * pz1;
        itsMap[PY] = euclid.M(0, 1) * px1 + euclid.M(1, 1) * py1 + euclid.M(2, 1) * pz1;
        Linear pz2 = euclid.M(0, 2) * px1 + euclid.M(1, 2) * py1 + euclid.M(2, 2) * pz1;

        Linear x = itsMap[X] - euclid.getX();
        Linear y = itsMap[Y] - euclid.getY();
        Linear x2 =
            euclid.M(0, 0) * x + euclid.M(1, 0) * y - euclid.M(2, 0) * euclid.getZ();
        Linear y2 =
            euclid.M(0, 1) * x + euclid.M(1, 1) * y - euclid.M(2, 1) * euclid.getZ();
        Linear s2 =
            euclid.M(0, 2) * x + euclid.M(1, 2) * y - euclid.M(2, 2) * euclid.getZ();
        Linear sByPz = s2 / pz2;

        double kin = itsReference.getM() / itsReference.getP();
        Linear E = sqrt(pt * pt + kin * kin);
        double refTime = refLength / itsReference.getBeta();
        itsMap[X] = x2 - sByPz * itsMap[PX];
        itsMap[Y] = y2 - sByPz * itsMap[PY];
        itsMap[T] += pt * (refTime / E  + sByPz);
    }
}
