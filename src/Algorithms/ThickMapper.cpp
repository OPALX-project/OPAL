// ------------------------------------------------------------------------
// $RCSfile: ThickMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickMapper
//   The visitor class for building a map of given order for a beamline
//   using a finite-length lenses for all elements.
//   Multipole-like elements are done by expanding the Lie series.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/ThickMapper.h"

#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"

#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FVps.h"

#include "Physics/Physics.h"

#include "Utilities/NumToStr.h"

class Beamline;
class PartData;
using Physics::c;

#define PSdim 6
typedef FVector<double, PSdim> Vector;
typedef FMatrix<double, PSdim, PSdim> Matrix;
typedef FTps<double, PSdim> Series;
typedef FVps<double, PSdim> Map, VSeries;
typedef FMatrix<FTps<double, PSdim>, PSdim, PSdim> MxSeries;

namespace {
    Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds,
                           int nx = 20);
    Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds,
                        int nx = 20, int cx = 4);
};

// Class ThickMapper
// ------------------------------------------------------------------------

ThickMapper::ThickMapper(const Beamline &beamline,
                         const PartData &reference,
                         bool backBeam, bool backTrack):
    Mapper(beamline, reference, backBeam, backTrack)
{}


ThickMapper::~ThickMapper()
{}


void ThickMapper::visitBeamBeam(const BeamBeam &) {
    // *** MISSING *** Map for beam-beam.
}

void ThickMapper::visitBeamStripping(const BeamStripping &) {
    // *** MISSING *** Map for beam Stripping.
}

void ThickMapper::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}


void ThickMapper::visitCorrector(const Corrector &corr) {
    double length = flip_s * corr.getElementLength();

    // Drift through first half of length.
    if(length != 0.0) applyDrift(length / 2.0);

    // Apply kick.
    double scale =
        (flip_s * flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BDipoleField &field = corr.getField();
    itsMap[PX] -= field.getBy() * scale;
    itsMap[PY] += field.getBx() * scale;

    // Drift through second half of length.
    if(length != 0.0) applyDrift(length / 2.0);
}

void ThickMapper::visitDegrader(const Degrader &deg) {
    applyDrift(flip_s * deg.getElementLength());
}

void ThickMapper::visitDiagnostic(const Diagnostic &diag) {
    // The diagnostic has no effect on the map.
    applyDrift(flip_s * diag.getElementLength());
}


void ThickMapper::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void ThickMapper::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void ThickMapper::visitLambertson(const Lambertson &lamb) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void ThickMapper::visitMarker(const Marker &marker) {
    // Do nothing.
}


void ThickMapper::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void ThickMapper::visitMultipole(const Multipole &mult) {
    //std::cerr << "==> In ThickMapper::visitMultipole(const Multipole &mult)" << std::endl;
    // Geometry and Field
    double length = flip_s * mult.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = mult.getField();
    int order = field.order();

    if(length == 0.0) {
        // Thin multipole, field coefficients are integral(B*dl).
        scale *= flip_s;
        applyThinMultipole(field, scale);
    } else {
        // Finite-length multipole, field coefficients are B.
        // Apply entrance fringe field.
        applyEntranceFringe(0.0, 0.0, field, scale);

        // Remove closed orbit from map.
        Vector co = itsMap.constantTerm();
        itsMap.setMinOrder(1);

        // Construct the Hamiltonian H about that closed orbit.
        // (1) Define variables.
        static const Series px = Series::makeVariable(PX);
        static const Series py = Series::makeVariable(PY);
        static const Series pt = Series::makeVariable(PT) + 1.0; // 1 + \delta = p/p0
        static const Map ident;
        Map translate = co + ident;
        // (2) Construct kinematic terms.
        double kin = itsReference.getM() / itsReference.getP(); // = 1/(beta*gamma)
        Series E_trans = sqrt((pt * pt + kin * kin).substitute(translate), order) / itsReference.getBeta();
        Series pz_trans  = sqrt((pt * pt - px * px - py * py).substitute(translate), order);
        // (3) Build vector potential in straight reference frame.
        Series As = buildMultipoleVectorPotential(field) * scale;
        As.setTruncOrder(Series::EXACT);
        // (4) Build Hamiltonian.
        Series H_trans = E_trans - pz_trans + As.substitute(translate);

        // Propagate closed orbit.
        // Build J.grad(H).
        VSeries JgradH;
        for(int i = 0; i < PSdim; i += 2) {
            JgradH[ i ] =  H_trans.derivative(i + 1);
            JgradH[i+1] = -H_trans.derivative(i);
        }
        // Do integration to propagate closed orbit using implicitInt4().
        static const Vector zeroV;
        Vector new_co = co + implicitInt4(zeroV, JgradH, length, 0.5 * length);

        // Propagate map.
        H_trans.setMinOrder(2); // Can do this because we assume zin is on the closed orbit.
        itsMap = ExpMap((-length) * H_trans).substitute(itsMap);

        // Add new closed orbit to propagated map.
        itsMap += new_co;

        // Apply exit fringe field.
        applyExitFringe(0.0, 0.0, field, scale);
    }
    //std::cerr << "==> Leaving ThickMapper::visitMultipole(...)" << std::endl;
}

void ThickMapper::visitProbe(const Probe &prob) {
    // Do nothing.
}

void ThickMapper::visitRBend(const RBend &bend) {
    //std::cerr << "==> In ThickMapper::visitRBend(const RBend &bend)" << std::endl;
    // Geometry and Field.
    const RBendGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();
    int order = field.order();
    double beta_inv = 1.0 / itsReference.getBeta();

    // Compute slices and stepsize.
    int nSlices = (int) bend.getSlices();
    double stepsize = bend.getStepsize();
    if(stepsize != 0.0) {
        int nst = (int) ceil(length / stepsize);
        nSlices = std::max(nSlices, nst);
    }
    double d_ell = length / nSlices;

    if(length == 0.0) {
        // Thin RBend.
        double half_angle = flip_s * geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);
        applyTransform(rotat, 0.0); // Transform from in-plane to mid-plane.
        applyThinMultipole(field, scale); // Apply multipole kick.
        applyTransform(rotat, 0.0); // Transform from mid-plane to out-plane.
    } else {
        // Finite-length RBend.
        // Define variables, ...
        static const Series x  = Series::makeVariable(X);
        static const Series px = Series::makeVariable(PX);
        static const Series py = Series::makeVariable(PY);
        static const Series de = Series::makeVariable(PT);
        static const Series de2 = de * de;
        static const Series pxy2 = px * px + py * py;
        static const Map ident;
        //   kinematic terms, ...
        double kin = itsReference.getM() / itsReference.getP(); // mc^2/pc = 1/(beta*gamma)
        double kin2 = kin * kin;
        //   and the vector potential, actually -As in straight reference frame.
        Series As = buildMultipoleVectorPotential(field) * scale;
        As.setTruncOrder(Series::EXACT);

        if(back_path) {
            // Apply rotation global to local.
            applyTransform(Inverse(geometry.getExitPatch()));
            // Apply entrance fringe field.
            applyEntranceFringe(bend.getExitFaceRotation(), scale, field, scale);
        } else {
            // Apply rotation global to local.
            applyTransform(geometry.getEntrancePatch());
            // Apply exit fringe field.
            applyEntranceFringe(bend.getEntryFaceRotation(), scale, field, scale);
        }

        // Remove closed orbit from map.
        Vector co = itsMap.constantTerm();
        itsMap.setMinOrder(1);

        // Propagate closed orbit and map.
        for(int slice = 0; slice < nSlices; ++slice) {
            // Construct translated Hamiltonian about closed orbit.
            double p1 = 1.0 + co[PT];
            Series de12 = p1 * p1 + 2.0 * p1 * de + de2;
            Series H_trans = beta_inv * sqrt(de12 + kin2, order)
                             - sqrt(de12 - co[PX] * co[PX] - co[PY] * co[PY]
                                    - 2.0 * (co[PX] * px + co[PY] * py) - pxy2, order)
                             + As.substitute(co + ident);
            // Build J.grad(H).
            VSeries JgradH;
            for(int i = 0; i < PSdim; i += 2) {
                JgradH[ i ] =  H_trans.derivative(i + 1);
                JgradH[i+1] = -H_trans.derivative(i);
            }
            // Do integration to propagate closed orbit using implicitInt4().
            static const Vector zeroV;
            Vector new_co = co + implicitInt4(zeroV, JgradH, d_ell, 0.5 * d_ell);
            // Propagate map using mid-point of closed-orbit path.
            Vector mid_co = (co + new_co) / 2.0;
            p1 = 1.0 + mid_co[PT];
            de12 = p1 * p1 + 2.0 * p1 * de + de2;
            H_trans = beta_inv * sqrt(de12 + kin2, order)
                      - sqrt(de12 - mid_co[PX] * mid_co[PX] - mid_co[PY] * mid_co[PY]
                             - 2.0 * (mid_co[PX] * px + mid_co[PY] * py) - pxy2, order)
                      + As.substitute(mid_co + ident);
            H_trans.setMinOrder(2); // Can do this because we expand about the closed orbit.
            itsMap = ExpMap((-d_ell) * H_trans).substitute(itsMap);
            // Set closed-orbit to new value.
            co = new_co;
        }

        // Add new closed orbit to propagated map.
        itsMap += co;

        if(back_path) {
            // Apply exit fringe field.
            applyExitFringe(bend.getEntryFaceRotation(), scale, field, scale);
            // Transform from local to global.
            applyTransform(Inverse(geometry.getEntrancePatch()));
        } else {
            // Apply entrance fringe field.
            applyEntranceFringe(bend.getExitFaceRotation(), scale, field, scale);
            // Transform from local to global.
            applyTransform(geometry.getExitPatch());
        }
    }
}


void ThickMapper::visitRFCavity(const RFCavity &as) {
    double length = flip_s * as.getElementLength();

    // Drift through half length.
    if(length != 0.0) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    double freq = as.getFrequency();
    double kin = itsReference.getM() / itsReference.getP();
    double peak = flip_s * as.getAmplitude() / itsReference.getP();

    Series pt = itsMap[PT] + 1.0;
    Series speed = (c * pt) / sqrt(pt * pt + kin * kin);
    Series phase = as.getPhase() + freq * itsMap[T] / speed;
    itsMap[PT] += peak * sin(phase) / pt;

    // Drift through half length.
    if(length != 0.0) applyDrift(length / 2.0);
}


void ThickMapper::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Map for RF quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void ThickMapper::visitSBend(const SBend &bend) {
    //std::cerr << "==> In ThickMapper::visitSBend(const SBend &bend)" << std::endl;
    // Geometry and Field.
    const PlanarArcGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();
    int order = field.order();
    double beta_inv = 1.0 / itsReference.getBeta();

    // Compute slices and stepsize.
    int nSlices = (int) bend.getSlices();
    double stepsize = bend.getStepsize();
    if(stepsize != 0.0) {
        int nst = (int) ceil(length / stepsize);
        nSlices = std::max(nSlices, nst);
    }
    double d_ell = length / nSlices;

    if(length == 0.0) {
        // Thin SBend.
        double half_angle = geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);
        applyTransform(rotat, 0.0); // Transform from in-plane to mid-plane.
        applyThinSBend(field, scale, 0.0); // Apply multipole kick.
        applyTransform(rotat, 0.0); // Transform from mid-plane to out-plane.
    } else {
        // Finite-length SBend.
        // Define variables, ...
        static const Series x  = Series::makeVariable(X);
        static const Series px = Series::makeVariable(PX);
        static const Series py = Series::makeVariable(PY);
        static const Series de = Series::makeVariable(PT);
        static const Series de2 = de * de;
        static const Series pxy2 = px * px + py * py;
        static const Map ident;
        //   kinematic terms, ...
        double h = geometry.getCurvature();
        double kin = itsReference.getM() / itsReference.getP(); // mc^2/pc = 1/(beta*gamma)
        double kin2 = kin * kin;
        Series hx1 = (1.0 + h * x);
        //   and the vector potential, actually -(1+h*x)*As in reference frame of curvature h.
        Series As = buildSBendVectorPotential(field, h) * scale;
        As.setTruncOrder(Series::EXACT);

        // Apply entrance fringe field.
        applyEntranceFringe(bend.getEntryFaceRotation(),
                            bend.getEntryFaceCurvature(), field, scale);

        // Remove closed orbit from map.
        Vector co = itsMap.constantTerm();
        itsMap.setMinOrder(1);

        // Propagate closed orbit and map.
        for(int slice = 0; slice < nSlices; ++slice) {
            // Construct translated Hamiltonian about closed orbit.
            double p1 = 1.0 + co[PT];
            Series de12 = p1 * p1 + 2.0 * p1 * de + de2;
            Series H_trans = beta_inv * sqrt(de12 + kin2, order)
                             - (hx1 + h * co[X])
                             * sqrt(de12 - co[PX] * co[PX] - co[PY] * co[PY]
                                    - 2.0 * (co[PX] * px + co[PY] * py) - pxy2, order)
                             + As.substitute(co + ident);
            // Build J.grad(H).
            VSeries JgradH;
            for(int i = 0; i < PSdim; i += 2) {
                JgradH[ i ] =  H_trans.derivative(i + 1);
                JgradH[i+1] = -H_trans.derivative(i);
            }
            // Do integration to propagate closed orbit using implicitInt4().
            static const Vector zeroV;
            Vector new_co = co + implicitInt4(zeroV, JgradH, d_ell, 0.5 * d_ell);
            // Propagate map using mid-point of closed-orbit path.
            Vector mid_co = (co + new_co) / 2.0;
            p1 = 1.0 + mid_co[PT];
            de12 = p1 * p1 + 2.0 * p1 * de + de2;
            H_trans = beta_inv * sqrt(de12 + kin2, order)
                      - (hx1 + h * mid_co[X])
                      * sqrt(de12 - mid_co[PX] * mid_co[PX] - mid_co[PY] * mid_co[PY]
                             - 2.0 * (mid_co[PX] * px + mid_co[PY] * py) - pxy2, order)
                      + As.substitute(mid_co + ident);
            H_trans.setMinOrder(2); // Can do this because we expand about the closed orbit.
            itsMap = ExpMap((-d_ell) * H_trans).substitute(itsMap);
            // Set closed-orbit to new value.
            co = new_co;
        }

        // Add new closed orbit to propagated map.
        itsMap += co;

        // Apply exit fringe field.
        applyExitFringe(bend.getExitFaceRotation(),
                        bend.getExitFaceCurvature(), field, scale);
    }
    //std::cerr << "==> Leaving ThickMapper::visitSBend(...)" << std::endl;
}


void ThickMapper::visitSeparator(const Separator &sep) {
    double length = flip_s * sep.getElementLength();

    if(length != 0.0) {
        // Drift through first half of length.
        applyDrift(length / 2.0);

        // electrostatic kick.
        double scale = (length * itsReference.getQ()) / itsReference.getP();
        double Ex = scale * sep.getEx();
        double Ey = scale * sep.getEy();
        Series pt = 1.0 + itsMap[PT];
        itsMap[PX] += Ex / pt;
        itsMap[PY] += Ey / pt;

        // Drift through second half of length.
        applyDrift(length / 2.0);
    }
}


void ThickMapper::visitSeptum(const Septum &sept) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void ThickMapper::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length != 0.0) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double kin = itsReference.getM() / itsReference.getP();
            double refTime = length / itsReference.getBeta();

            Series pt = itsMap[PT] + 1.0;
            Series px = itsMap[PX] + ks * itsMap[Y];
            Series py = itsMap[PY] - ks * itsMap[X];
            Series pz = sqrt(pt * pt - px * px - py * py);
            Series E = sqrt(pt * pt + kin * kin);

            Series k = ks / pz;
            Series C = cos(k * length);
            Series S = sin(k * length);

            Series xt  = C * itsMap[X]  + S * itsMap[Y];
            Series yt  = C * itsMap[Y]  - S * itsMap[X];
            Series pxt = C * itsMap[PX] + S * itsMap[PY];
            Series pyt = C * itsMap[PY] - S * itsMap[PX];

            itsMap[X]  = C * xt  + (S / k) * pxt;
            itsMap[Y]  = C * yt  + (S / k) * pyt;
            itsMap[PX] = C * pxt - (S * k) * xt;
            itsMap[PY] = C * pyt - (S * k) * yt;
            itsMap[T] += pt * (refTime / E - length / pz);
        } else {
            applyDrift(length);
        }
    }
}


void ThickMapper::visitParallelPlate(const ParallelPlate &pplate) {
    // Do nothing.
}


void ThickMapper::applyDrift(double length) {
    if(length != 0.0) {
        // HACK: fix a maximum order for computations.
        int order = 6;

        // Remove closed orbit from map.
        Vector co = itsMap.constantTerm();
        itsMap.setMinOrder(1);

        // Construct the Hamiltonian H about that closed orbit.
        // (1) Define variables.
        static const Series px = Series::makeVariable(PX);
        static const Series py = Series::makeVariable(PY);
        static const Series pt = Series::makeVariable(PT) + 1.0;
        static const Series x  = Series::makeVariable(X);
        static const Map ident;
        Map translate = co + ident;
        // (2) Construct kinematic terms.
        double kin = itsReference.getM() / itsReference.getP(); // 1/(beta*gamma)
        Series E_trans = sqrt((pt * pt + kin * kin).substitute(translate), order)
                         / itsReference.getBeta();
        Series pz_trans = sqrt((pt * pt - px * px - py * py).substitute(translate), order);
        // (3) Build Hamiltonian.
        Series H_trans = E_trans - pz_trans;

        // Propagate closed orbit.
        // Build J.grad(H).
        VSeries JgradH;
        for(int i = 0; i < PSdim; i += 2) {
            JgradH[ i ] =  H_trans.derivative(i + 1);
            JgradH[i+1] = -H_trans.derivative(i);
        }
        // Do integration to propagate closed orbit using implicitInt4().
        static const Vector zeroV;
        Vector new_co = co + implicitInt4(zeroV, JgradH, length, 0.5 * length);

        // Propagate map.
        H_trans.setMinOrder(2); // Can do this because we assume zin is on the closed orbit.
        itsMap = ExpMap((-length) * H_trans).substitute(itsMap);

        // Add new closed orbit to propagated map.
        itsMap += new_co;
    }
}


void ThickMapper::applyEntranceFringe(double angle, double curve,
                                      const BMultipoleField &field, double scale) {
    // *** MISSING *** Higher order terms for entrance fringe.
    double ca = cos(angle);
    double sa = sin(angle);
    double ta = tan(angle);

    int order = field.order();
    static const Series x = Series::makeVariable(X);
    Series by = field.normal(order);
    for(int i = order; --i >= 1;) by = by * x + field.normal(i);
    by *= scale;

    const Series p = itsMap[PT] + 1.0;
    // rotate to magnet face
    Series ps = sqrt(p * p - itsMap[PX] * itsMap[PX] - itsMap[PY] * itsMap[PY], order);
    itsMap[X] = itsMap[X] / (ca * (1.0 - ta * itsMap[PX] / ps));
    itsMap[PX] = ca * itsMap[PX] + sa * ps;
    Series ellovpp = sa * itsMap[X] / ps;
    itsMap[Y] += ellovpp * itsMap[PY];
    itsMap[T] -= ellovpp * p;
    // fringe
    Series psy = sqrt(p * p - itsMap[PX] * itsMap[PX], order);
    itsMap[PY] -= by.substitute(itsMap) * itsMap[Y] * itsMap[PX] / psy;
    // rotate from magnet face
    ps = sqrt(p * p - itsMap[PX] * itsMap[PX] - itsMap[PY] * itsMap[PY], order);
    itsMap[X] = itsMap[X] / (ca * (1.0 + ta * itsMap[PX] / ps));
    itsMap[PX] = ca * itsMap[PX] - sa * ps;
    ellovpp = sa * itsMap[X] / ps;
    itsMap[Y] -= ellovpp * itsMap[PY];
    itsMap[T] += ellovpp * p;
    // edge effect
    itsMap[PX] += by.substitute(itsMap) * ta * itsMap[X];
}


void ThickMapper::applyExitFringe(double angle, double curve,
                                  const BMultipoleField &field, double scale) {
    // *** MISSING *** Higher order terms for exit fringe.
    double ca = cos(angle);
    double sa = sin(angle);
    double ta = tan(angle);

    int order = field.order();
    static const Series x = Series::makeVariable(X);
    Series by = field.normal(order);
    for(int i = order; --i >= 1;) by = by * x + field.normal(i);
    by *= scale;

    const Series p = itsMap[PT] + 1.0;
    // edge effect
    itsMap[PX] += by.substitute(itsMap) * ta * itsMap[X];
    // rotate to magnet face
    Series ps = sqrt(p * p - itsMap[PX] * itsMap[PX] - itsMap[PY] * itsMap[PY], order);
    itsMap[X] = itsMap[X] / (ca * (1.0 + ta * itsMap[PX] / ps));
    itsMap[PX] = ca * itsMap[PX] - sa * ps;
    Series ellovpp = sa * itsMap[X] / ps;
    itsMap[Y] -= ellovpp * itsMap[PY];
    itsMap[T] += ellovpp * p;
    // fringe
    Series psy = sqrt(p * p - itsMap[PX] * itsMap[PX], order);
    itsMap[PY] += by.substitute(itsMap) * itsMap[Y] * itsMap[PX] / psy;
    // rotate from magnet face
    ps = sqrt(p * p - itsMap[PX] * itsMap[PX] - itsMap[PY] * itsMap[PY], order);
    itsMap[X] = itsMap[X] / (ca * (1.0 - ta * itsMap[PX] / ps));
    itsMap[PX] = ca * itsMap[PX] + sa * ps;
    ellovpp = sa * itsMap[X] / ps;
    itsMap[Y] += ellovpp * itsMap[PY];
    itsMap[T] -= ellovpp * p;
}

namespace {
    Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx) {
        //std::cerr << "==> In implicitInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
        // Default: nx = 20, cx = 4

        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
        // version of implicitInt2 to obtain zf accurate through fourth-order in the
        // step-size ds.  When f derives from a Hamiltonian---i.e., f = J.grad(H)---
        // then this routine performs symplectic integration.  The optional arguments
        // nx and cx have the same meaning as in implicitInt2().

        // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
        static const double yt = pow(2., 1 / 3.);
        static const double ya = 1 / (2. - yt);
        static const double yb = -yt * ya;

        // Build matrix grad(f).
        MxSeries gradf;
        for(int i = 0; i < PSdim; ++i)
            for(int j = 0; j < PSdim; ++j)
                gradf[i][j] = f[i].derivative(j);

        // Initialize accumulated length, current step-size, and number of cuts.
        double as = std::abs(s), st = 0., dsc = std::abs(ds);
        if(s < 0.) dsc = -dsc;
        int ci = 0;

        // Integrate each step.
        Vector zf = zin;
        while(std::abs(st) < as) {
            Vector zt;
            bool ok = true;
            try {
                if(std::abs(st + dsc) > as) dsc = s - st;
                zt = ::implicitIntStep(zf, f, gradf, ya * dsc, nx);
                zt = ::implicitIntStep(zt, f, gradf, yb * dsc, nx);
                zt = ::implicitIntStep(zt, f, gradf, ya * dsc, nx);
            } catch(ConvergenceError &cnverr) {
                if(++ci > cx) {
                    std::string msg = "Convergence not achieved within " + NumToStr<int>(cx) + " cuts of step-size!";
                    throw ConvergenceError("ThickMapper::implicitInt4()", msg);
                }
                //std::cerr << "  cutting step size in half" << std::endl;
                dsc *= 0.5;
                ok = false;
            }
            if(ok) {zf = zt; st += dsc;}
        }

        //std::cerr << "==> Leaving implicitInt4(zin,f,s,ds,nx,cx)" << std::endl;
        return zf;
    }

    Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds, int nx) {
        //std::cerr << "==> In implicitIntStep(zin,f,gradf,ds,nx) ..." << std::endl;
        //std::cerr << "  ds = " << ds << std::endl;
        //std::cerr << " zin =\n" << zin << std::endl;
        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a single step of size ds, using Newton's method to solve the
        // implicit equation zf = zin + ds*f((zin+zf)/2).  For reasons of efficiency,
        // its arguments include the matrix gradf = grad(f).  The (optional) argument
        // nx limits the number of Newton iterations.  This routine returns a result
        // zf accurate through second-order in the step-size ds.  When f derives from
        // a Hamiltonian---i.e., f=J.grad(H)---then this routine performs symplectic
        // integration.

        // Set up flags, etc., for convergence (bounce) test.
        FVector<bool, PSdim> bounce(false);
        Vector dz, dz_old;
        int bcount = 0;
        static const double thresh = 1.e-8;

        // Use second-order Runge-Kutta integration to determine a good initial guess.
        double ds2 = 0.5 * ds;
        Vector z = f.constantTerm(zin);
        z = zin + ds2 * (z + f.constantTerm(zin + ds * z));

        // Newton iterations:
        //   z :-> [I-ds/2.grad(f)]^{-1}.[zin+ds.f((zin+z)/2)-ds/2.grad(f).z]
        Vector zf;
        int ni = 0;
        while(bcount < PSdim) {
            if(ni == nx) {
                std::string msg = "Convergence not achieved within " + NumToStr<int>(nx) + " iterations!";
                throw ConvergenceError("ThickMapper::implicitIntStep()", msg);
            }

            // Build gf = -ds/2.grad(f)[(zin+z)/2] and idgf_inv = [I-ds/2.grad(f)]^{-1}[(zin+z)/2].
            Vector zt = 0.5 * (zin + z);
            Matrix gf, idgf, idgf_inv;
            for(int i = 0; i < PSdim; ++i)
                for(int j = 0; j < PSdim; ++j)
                    gf[i][j] = -ds2 * gradf[i][j].evaluate(zt);
            idgf = gf;
            for(int i = 0; i < PSdim; ++i) idgf[i][i] += 1.;
            FLUMatrix<double, PSdim> lu(idgf);
            idgf_inv = lu.inverse();

            // Execute Newton step.
            zf = idgf_inv * (zin + ds * f.constantTerm(zt) + gf * z);

            //std::cerr << " -(ds/2)grad(f) =\n" << gf << std::endl;
            //std::cerr << " f =\n" << f.constantTerm(zt) << std::endl;
            //std::cerr << "zk =\n" << zf << std::endl;

            // Test for convergence ("bounce" test).
            dz_old = dz;
            dz = zf - z;
            if(ni) { // (we need at least two iterations before testing makes sense)
                for(int i = 0; i < PSdim; ++i) {
                    if(!bounce[i] && (dz[i] == 0. || (std::abs(dz[i]) < thresh && std::abs(dz[i]) >= std::abs(dz_old[i]))))
                        {bounce[i] = true; ++bcount;}
                }
            }
            z = zf;
            ++ni;
        }

        //std::cerr << "  zf =\n" << zf << std::endl;
        //std::cerr << "==> Leaving implicitIntStep(zin,f,gradf,ds,nx)" << std::endl;
        return zf;
    }
};