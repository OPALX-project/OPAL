// ------------------------------------------------------------------------
// $RCSfile: ThinMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThinMapper
//   The visitor class for building a map for a beamline
//   using a thin-lens approximation for all elements.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/ThinMapper.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Undulator.h"
#include "AbsBeamline/Degrader.h"
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
#include "AbsBeamline/CyclotronValley.h"

#include "Algorithms/MapIntegrator.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "Physics/Physics.h"

using Physics::c;

typedef FTps<double, 6> Series;


// Class ThinMapper
// ------------------------------------------------------------------------

ThinMapper::ThinMapper(const Beamline &beamline, const PartData &ref,
                       bool backBeam, bool backTrack):
    Mapper(beamline, ref, backBeam, backTrack)
{}


ThinMapper::~ThinMapper()
{}


void ThinMapper::visitBeamBeam(const BeamBeam &bb) {
    // *** MISSING *** Map algorithm on BeamBeam
}


void ThinMapper::visitDegrader(const Degrader &deg) {
    applyDrift(flip_s * deg.getElementLength());
}

void ThinMapper::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}


void ThinMapper::visitCorrector(const Corrector &corr) {
    // Drift through first half of length.
    double length = flip_s * corr.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply kick.
    double scale = (flip_s * flip_B * itsReference.getQ() * c) /
                   itsReference.getP();
    const BDipoleField &field = corr.getField();
    itsMap[PX] -= field.getBy() * scale;
    itsMap[PY] += field.getBx() * scale;

    // Drift through second half of length.
    if(length) applyDrift(length / 2.0);
}


void ThinMapper::visitDiagnostic(const Diagnostic &diag) {
    // The diagnostic has no effect on the map.
    applyDrift(flip_s * diag.getElementLength());
}


void ThinMapper::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void ThinMapper::visitUndulator(const Undulator &u) {
     INFOMSG("Undulator not implemented");
}

void ThinMapper::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void ThinMapper::visitLambertson(const Lambertson &lamb) {
    // Assume that the reference orbit is in the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void ThinMapper::visitMarker(const Marker &) {
    // Do nothing.
}


void ThinMapper::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void ThinMapper::visitMultipole(const Multipole &mult) {
    double length = flip_s * mult.getElementLength();
    const BMultipoleField &field = mult.getField();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();

    if(length) {
        // Drift through first half of the length.
        applyDrift(length / 2.0);

        // Apply thin multipole kick, field is per unit length.
        scale *= length;
        applyThinMultipole(field, scale);

        // Drift through second half of the length.
        applyDrift(length / 2.0);
    } else {
        // Thin multipole, field is integral(K*dl).
        scale *= flip_s;
        applyThinMultipole(field, scale);
    }
}

void ThinMapper::visitProbe(const Probe &) {
    // Do nothing.
}


void ThinMapper::visitRBend(const RBend &bend) {
    // Geometry.
    const RBendGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double angle = flip_s * geometry.getBendAngle();

    // Magnetic field.
    const BMultipoleField &field = bend.getField();

    // Drift to mid-plane.
    applyDrift(length / 2.0);

    // Apply multipole kick and linear approximation to geometric bend.
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    if(length) scale *= length;
    int order = field.order();

    if(order > 0) {
        Series x = itsMap[X];
        Series y = itsMap[Y];
        Series kx = + field.normal(order);
        Series ky = - field.skew(order);

        while(--order > 0) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        itsMap[PX] -= kx * scale - angle * (1.0 + itsMap[PT]);
        itsMap[PY] += ky * scale;
        itsMap[T]  -= angle * x;
    }

    // Drift to out-plane.
    applyDrift(length / 2.0);
}


void ThinMapper::visitRFCavity(const RFCavity &as) {
    // Drift through half length.
    double length = flip_s * as.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    double kin = itsReference.getM() / itsReference.getP();
    double freq = as.getFrequency();
    double peak = flip_s * as.getAmplitude() / itsReference.getP();

    Series pt = itsMap[PT] + 1.0;
    Series speed = (c * pt) / sqrt(pt * pt + kin * kin);
    Series phase = as.getPhase() + freq * itsMap[T] / speed;
    itsMap[PT] += peak * sin(phase) / pt;

    // Drift through half length.
    if(length) applyDrift(length / 2.0);
}


void ThinMapper::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Map algorithm on RF Quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void ThinMapper::visitSBend(const SBend &bend) {
    const PlanarArcGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double angle = flip_s * geometry.getBendAngle();

    // Magnetic field.
    const BMultipoleField &field = bend.getField();

    // Drift to mid-plane.
    applyDrift(length / 2.0);

    // Apply multipole kick and linear approximation to geometric bend.
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    if(length) scale *= length;
    int order = field.order();

    if(order > 0) {
        Series x = itsMap[X];
        Series y = itsMap[Y];
        Series kx = + field.normal(order);
        Series ky = - field.skew(order);

        while(--order > 0) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        itsMap[PX] -= kx * scale - angle * (1.0 + itsMap[PT]);
        itsMap[PY] += ky * scale;
        itsMap[T]  -= angle * x;
    }

    // Drift to out-plane.
    applyDrift(length / 2.0);
}


void ThinMapper::visitSeparator(const Separator &sep) {
    // Drift through first half of length.
    double length = flip_s * sep.getElementLength();

    if(length) {
        applyDrift(length / 2.0);

        double scale = (length * itsReference.getQ()) / itsReference.getP();
        double Ex = scale * sep.getEx();
        double Ey = scale * sep.getEy();
        Series pt = 1.0 + itsMap[PT];
        itsMap[PX] += Ex / pt;
        itsMap[PY] += Ey / pt;

        applyDrift(length / 2.0);
    }
}


void ThinMapper::visitSeptum(const Septum &sept) {
    // Assume that the reference orbit is in the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void ThinMapper::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            Series pt = itsMap[PT] + 1.0;
            Series px = itsMap[PX] + ks * itsMap[Y];
            Series py = itsMap[PY] - ks * itsMap[X];
            Series pz = sqrt(pt * pt - px * px - py * py);
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

            double kin = itsReference.getM() / itsReference.getP();
            double ref = kin * kin;
            itsMap[T] += length * (pt * ref - (px * px + py * py + 3.0 * pt * pt * ref) / 2.0);
        } else {
            applyDrift(length);
        }
    }
}


void ThinMapper::visitParallelPlate(const ParallelPlate &) {
    // Do nothing.
}

void ThinMapper::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing.
}

void ThinMapper::applyDrift(double length) {
    double kin  = itsReference.getM() / itsReference.getP();
    double ref  = kin * kin;
    Series px = itsMap[PX];
    Series py = itsMap[PY];
    Series pt = itsMap[PT];
    Series lByPz = length / (1.0 + pt);
    itsMap[X] += px * lByPz;
    itsMap[Y] += py * lByPz;
    itsMap[T] += length * (pt * ref - (px * px + py * py + 3.0 * pt * pt * ref) / 2.0);
}
