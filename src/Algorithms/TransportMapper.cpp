// ------------------------------------------------------------------------
// $RCSfile: TransportMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.4 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TransportMapper
//   The visitor class for building a transport map for a beamline using
//   transport maps for all elements.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 08:16:30 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Algorithms/TransportMapper.h"

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
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
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/TransportFun.h"
#include "FixedAlgebra/TransportMath.h"
#include "Physics/Physics.h"
#include <cmath>

using Physics::c;

typedef FTps<double, 2> Series2;
typedef TransportFun<double, 6> TptFun;


// Class TransportMapper
// ------------------------------------------------------------------------

TransportMapper::TransportMapper(const Beamline &beamline,
                                 const PartData &reference,
                                 bool revBeam, bool revTrack):
    AbstractMapper(beamline, reference, revBeam, revTrack)
{}


TransportMapper::~TransportMapper()
{}


void TransportMapper::getMap(LinearMap<double, 6> &map) const {
    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j <= 6; ++j) {
            map[i][j] = itsMap[i][j];
        }
    }
}


void TransportMapper::getMap(TransportMap<double, 6> &map) const {
    map = itsMap;
}


void TransportMapper::getMap(FVps<double, 6> &map) const {
    map = FVps<double, 6>(itsMap);
}


void TransportMapper::setMap(const LinearMap<double, 6> &map) {
    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j <= 6; ++j) {
            itsMap[i][j] = map[i][j];
        }

        for(int j = 7; j < 28; ++j) {
            itsMap[i][j] = 0.0;
        }
    }
}


void TransportMapper::setMap(const TransportMap<double, 6> &map) {
    itsMap = map;
}


void TransportMapper::setMap(const FVps<double, 6> &map) {
    itsMap = TransportMap<double, 6>(map);
}


void TransportMapper::visitBeamBeam(const BeamBeam &) {
    // *** MISSING *** Map for beam-beam.
}

void TransportMapper::visitBeamStripping(const BeamStripping &) {
    // *** MISSING *** Map for beam Stripping.
}

void TransportMapper::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}


void TransportMapper::visitComponent(const Component &comp) {
    FVps<double, 6> map(itsMap);
    comp.trackMap(map, itsReference, back_beam, back_track);
    itsMap = TransportMap<double, 6>(map);
}


void TransportMapper::visitCorrector(const Corrector &corr) {
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

void TransportMapper::visitDegrader(const Degrader &deg) {
    applyDrift(flip_s * deg.getElementLength());
}

void TransportMapper::visitDiagnostic(const Diagnostic &diag) {
    applyDrift(flip_s * diag.getElementLength());
}


void TransportMapper::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void TransportMapper::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void TransportMapper::visitLambertson(const Lambertson &lamb) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void TransportMapper::visitMarker(const Marker &marker) {
    // Do nothing.
}


void TransportMapper::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void TransportMapper::visitMultipole(const Multipole &multipole) {
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


void TransportMapper::visitPatch(const Patch &patch) {
    Euclid3D transform = patch.getPatch();
    if(back_path) transform = Inverse(transform);
    applyTransform(transform, 0.0);
}
void TransportMapper::visitProbe(const Probe &prob) {
    // Do nothing.
}


void TransportMapper::visitRBend(const RBend &bend) {
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


void TransportMapper::visitRFCavity(const RFCavity &as) {
    // Drift through half length.
    double length = flip_s * as.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    TptFun time = itsMap[T] / (c * itsReference.getBeta());
    double peak = flip_s * as.getAmplitude() / itsReference.getP();
    TptFun phase = as.getPhase() + as.getFrequency() * time;
    itsMap[PT] += peak * sin(phase);

    // Drift through half length.
    if(length) applyDrift(length / 2.0);
}


void TransportMapper::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Map for RF Quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void TransportMapper::visitSBend(const SBend &bend) {
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


void TransportMapper::visitSeparator(const Separator &sep) {
    // Drift through first half of length.
    double length = flip_s * sep.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Electrostatic kick.
    double scale = (length * itsReference.getQ()) / itsReference.getP();
    double Ex = scale * sep.getEx();
    double Ey = scale * sep.getEy();
    TptFun pt = 1.0 + itsMap[PT];
    itsMap[PX] += Ex / pt;
    itsMap[PY] += Ey / pt;

    if(length) applyDrift(length / 2.0);
}


void TransportMapper::visitSeptum(const Septum &sept) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void TransportMapper::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double kin = itsReference.getM() / itsReference.getP();
            double refTime = length / itsReference.getBeta();
            TptFun pt = itsMap[PT] + 1.0;
            TptFun px = itsMap[PX] + ks * itsMap[Y];
            TptFun py = itsMap[PY] - ks * itsMap[X];
            TptFun pz = sqrt(pt * pt - px * px - py * py);
            TptFun E = sqrt(pt * pt + kin * kin);
            TptFun k = ks / pz;

            TptFun C = cos(k * length);
            TptFun S = sin(k * length);

            TptFun xt  = C * itsMap[X]  + S * itsMap[Y];
            TptFun yt  = C * itsMap[Y]  - S * itsMap[X];
            TptFun pxt = C * itsMap[PX] + S * itsMap[PY];
            TptFun pyt = C * itsMap[PY] - S * itsMap[PX];

            itsMap[X]  = C * xt  + (S / ks) * pxt;
            itsMap[Y]  = C * yt  + (S / ks) * pyt;
            itsMap[PX] = C * pxt - (S * ks) * xt;
            itsMap[PY] = C * pyt - (S * ks) * yt;
            itsMap[T] += pt * (refTime / E - length / pz);
        } else {
            applyDrift(length);
        }
    }
}


void TransportMapper::visitParallelPlate(const ParallelPlate &pplate) {
    // Do nothing.
}

void TransportMapper::visitAlignWrapper(const AlignWrapper &wrap) {
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


void TransportMapper::visitMapIntegrator(const MapIntegrator &i) {
    FVps<double, 6> map(itsMap);
    i.trackMap(map, itsReference, back_beam, back_track);
    itsMap = TransportMap<double, 6>(map);
}


// Helper function for finding first-order coefficients.
static void makeFocus
(double k, double L, double &c, double &s, double &d, double &f) {
    double t = k * L * L;
    if(std::abs(t) < 1.0e-4) {
        c = 1.0 - t / 2.0;
        s = L * (1.0 - t / 6.0);
        d = L * L * (0.5 - t / 24.0);
        f = L * L * L * (1.0 / 6.0 - t / 120.0);
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


void TransportMapper::applyDrift(double length) {
    double kin = itsReference.getM() / itsReference.getP();
    double refTime = length / itsReference.getBeta();

    TptFun px = itsMap[PX];
    TptFun py = itsMap[PY];
    TptFun pt = itsMap[PT] + 1.0;
    TptFun pz = sqrt(pt * pt - px * px - py * py);
    TptFun E = sqrt(pt * pt + kin * kin);

    itsMap[X] += length * px / pz;
    itsMap[Y] += length * py / pz;
    itsMap[T] += pt * (refTime / E - length / pz);
}


void TransportMapper::applyEntranceFringe(double angle,
        const BMultipoleField &field,
        double scale) {
    // ***** MISSING: second-order fringe at entrance
    double hx = scale * field.normal(1);
    double ex = hx * tan(angle);
    double ey = hx * tan(angle + itsMap[PX][0]);
    itsMap[PX] += ex * itsMap[X];
    itsMap[PY] -= ey * itsMap[Y];
}


void TransportMapper::applyExitFringe(double angle,
                                      const BMultipoleField &field,
                                      double scale) {
    // ***** MISSING: second-order fringe at exit
    double hx = scale * field.normal(1);
    double ex = hx * tan(angle);
    double ey = hx * tan(angle - itsMap[PX][0]);
    itsMap[PX] += ex * itsMap[X];
    itsMap[PY] -= ey * itsMap[Y];
}


void TransportMapper::applyTransportMap
(double length, double refLength, double h,
 const Series2 &Fx, const Series2 &Fy) {
    // ***** MISSING: Transport map for thick element.
    // Extract the phase coordinates.
    TptFun x  = itsMap[X];
    TptFun px = itsMap[PX];
    TptFun y  = itsMap[Y];
    TptFun py = itsMap[PY];
    TptFun pt = itsMap[PT];

    // (x0, y0) is the position of the entry point,
    // (x,  y)  is the position relative to the entry point.
    double x0 = x[0];
    double y0 = y[0];
    x.setCoefficient(0, 0.0);
    y.setCoefficient(0, 0.0);

    // Extract coefficients for equations of motion.
    // Indexing: 0 = constant, 1 = X, 2 = Y.
    double kx = Fx[1];
    double ks = (Fx[2] - Fy[1]) / 2.0;
    double ky = - Fy[2];
    double hx = h - Fx[0];
    double hy = Fy[0];

    // Longitudinal data.
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
        itsMap[X]  = cx * x + sx * px + dx * hx + h * dx * pt + x0;
        itsMap[PX] = wx * x + cx * px + sx * hx + h * sx * pt;
        itsMap[Y]  = cy * y + sy * py + dy * hy + y0;
        itsMap[PY] = wy * y + cy * py + sy * hy;
        itsMap[T] += refTime * pt + h * (sx * x + dx * px + fx * hx + length * x0);
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
        TptFun u  = c2 * x  - s2 * y;
        TptFun v  = c2 * y  + s2 * x;
        TptFun pu = c2 * px - s2 * py;
        TptFun pv = c2 * py + s2 * px;
        double hu = c2 * (h + hx) - s2 * hy;
        double hv = c2 * hy + s2 * (h + hx);
        double bu = c2 * h;
        double bv = s2 * h;

        // Advance through field.
        itsMap[X] = x0 +
                    ((cu + cv) * x + (cu - cv) * u + (su + sv) * px + (su - sv) * pu +
                     (du + dv) * hx + (du - dv) * hu + ((du + dv) * h + (du - dv) * bu) * pt) / 2.0;
        itsMap[PX] =
            ((wu + wv) * x + (wu - wv) * u + (cu + cv) * px + (cu - cv) * pu +
             (su + sv) * hx + (su - sv) * hu + ((su + sv) * h + (su - sv) * bu) * pt) / 2.0;
        itsMap[Y] = y0 +
                    ((cu + cv) * y - (cu - cv) * v + (su + sv) * py - (su - sv) * pv +
                     (du + dv) * hy - (du - dv) * hv - (du - dv) * bv * pt) / 2.0;
        itsMap[PY] =
            ((wu + wv) * y - (wu - wv) * v + (cu + cv) * py - (cu - cv) * pv +
             (su + sv) * hy - (su - sv) * hv - (su - sv) * bv * pt) / 2.0;
        itsMap[T] += refTime * pt;
    }
}


void TransportMapper::applyMultipoleBody(double length, double refLength,
        const BMultipoleField &field,
        double scale) {
    // Determine normalised field coefficients around actual orbit.
    // Fx and Fy are the normalised field coefficients,
    // expressed as transport functions in x and y.
    int order = field.order();
    Series2 Fx;
    Series2 Fy;
    Series2::setGlobalTruncOrder(2);

    if(order > 0) {
        static Series2 x = Series2::makeVariable(0);
        static Series2 y = Series2::makeVariable(1);
        x.setCoefficient(0, itsMap[X][0]);
        y.setCoefficient(0, itsMap[Y][0]);
        Fx = + field.normal(order);
        Fy = - field.skew(order);

        while(order > 1) {
            Series2 Fxt = x * Fx - y * Fy;
            Series2 Fyt = x * Fy + y * Fx;
            order--;
            Fx = Fxt + field.normal(order);
            Fy = Fyt - field.skew(order);
        };

        Fx *= scale;
        Fy *= scale;
    }

    applyTransportMap(length, refLength, 0.0, Fx, Fy);
}


void TransportMapper::
applySBendBody(double length, double refLength, double h,
               const BMultipoleField &field, double scale) {
    // Determine normalised field coefficients around actual orbit.
    // As is the vector potential times (1 + h*x),
    // expressed as a transport function in x and y.
    Series2 As = buildSBendVectorPotential(field, h) * scale;

    // Fx and Fy are the normalised field coefficients times (1 + h*x).
    Series2 Fx = + As.derivative(0);
    Series2 Fy = - As.derivative(1);
    applyTransportMap(length, refLength, h, Fx, Fy);
}


void TransportMapper::
applyThinMultipole(const BMultipoleField &field, double scale) {
    int order = field.order();

    if(order > 0) {
        TptFun x = itsMap[X];
        TptFun y = itsMap[Y];
        TptFun kx = + field.normal(order);
        TptFun ky = - field.skew(order);

        while(--order > 0) {
            TptFun kxt = x * kx - y * ky;
            TptFun kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        itsMap[PX] -= kx * scale;
        itsMap[PY] += ky * scale;
    }
}


void TransportMapper::
applyThinSBend(const BMultipoleField &field, double scale, double h) {
    Series2 As = buildSBendVectorPotential(field, h) * scale;

    // Fx and Fy are the normalised kicks,
    // expressed as transport functions in x and y.
    Series2 Fx = + As.derivative(0);
    Series2 Fy = - As.derivative(1);

    // These substitutions work because As depends on x and y only,
    // and not on px or py.
    itsMap[PX] -= Fx[0] + Fx[1] * itsMap[X] + Fx[2] * itsMap[Y];
    itsMap[PY] -= Fy[0] + Fy[1] * itsMap[X] + Fy[2] * itsMap[Y];
}


void TransportMapper::
applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        TptFun px1 = itsMap[PX];
        TptFun py1 = itsMap[PY];
        TptFun pt  = itsMap[PT] + 1.0;
        TptFun pz1 = sqrt(pt * pt - px1 * px1 - py1 * py1);

        itsMap[PX] = euclid.M(0, 0) * px1 + euclid.M(1, 0) * py1 + euclid.M(2, 0) * pz1;
        itsMap[PY] = euclid.M(0, 1) * px1 + euclid.M(1, 1) * py1 + euclid.M(2, 1) * pz1;
        TptFun pz2 = euclid.M(0, 2) * px1 + euclid.M(1, 2) * py1 + euclid.M(2, 2) * pz1;

        TptFun x = itsMap[X] - euclid.getX();
        TptFun y = itsMap[Y] - euclid.getY();
        TptFun x2 =
            euclid.M(0, 0) * x + euclid.M(1, 0) * y - euclid.M(2, 0) * euclid.getZ();
        TptFun y2 =
            euclid.M(0, 1) * x + euclid.M(1, 1) * y - euclid.M(2, 1) * euclid.getZ();
        TptFun s2 =
            euclid.M(0, 2) * x + euclid.M(1, 2) * y - euclid.M(2, 2) * euclid.getZ();
        TptFun sByPz = s2 / pz2;

        double kin = itsReference.getM() / itsReference.getP();
        TptFun E = sqrt(pt * pt + kin * kin);
        double refTime = refLength / itsReference.getBeta();
        itsMap[X] = x2 - sByPz * itsMap[PX];
        itsMap[Y] = y2 - sByPz * itsMap[PY];
        itsMap[T] += pt * (refTime / E  + sByPz);
    }
}


Series2 TransportMapper::
buildSBendVectorPotential(const BMultipoleField &field, double h) {
    int order = field.order();
    Series2::setGlobalTruncOrder(order + 1);
    Series2 As;

    if(order > 0) {
        static Series2 x = Series2::makeVariable(0);
        static Series2 y = Series2::makeVariable(1);
        x.setCoefficient(0, itsMap[X][0]);
        y.setCoefficient(0, itsMap[Y][0]);

        // Construct terms constant and transport in y.
        Series2 Ae = + field.normal(order); // Term even in y.
        Series2 Ao = - field.skew(order);   // Term odd  in y.
        int i = order;

        while(i > 1) {
            i--;
            Ae = Ae * x + field.normal(i);
            Ao = Ao * x - field.skew(i);
        };

        Ae = + (Ae * (1.0 + h * x)).integral(0);
        Ao = - (Ao * (1.0 + h * x));

        // Add terms up to maximum order.
        As = Ae + y * Ao;
        int k = 2;

        if(k <= order) {
            Series2 yp = y * y / 2.0;
            Series2 factor = h / (1.0 + h * x);

            while(true) {
                // Terms even in y.
                Ae = Ae.derivative(0);
                Ae = (factor * Ae - Ae.derivative(0)) * yp;
                As += Ae;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(0);
                Ao = (factor * Ao - Ao.derivative(0)) * yp;
                As += Ao;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    return As;
}