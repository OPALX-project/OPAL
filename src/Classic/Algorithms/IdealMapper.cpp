// ------------------------------------------------------------------------
// $RCSfile: IdealMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: IdealMapper
//   The visitor class for building a linear map for a beamline using
//   linear maps for all elements.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/IdealMapper.h"

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Patch.h"
#include "AbsBeamline/Separator.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "ComponentWrappers/CorrectorWrapper.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "ComponentWrappers/RBendWrapper.h"
#include "ComponentWrappers/SBendWrapper.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/LinearFun.h"
#include "FixedAlgebra/LinearMath.h"
#include "Physics/Physics.h"
#include <cmath>

using Physics::c;

typedef FTps<double, 2> Series2;
typedef LinearFun<double, 6> Linear;


// Class IdealMapper
// ------------------------------------------------------------------------

IdealMapper::IdealMapper(const Beamline &beamline, const PartData &reference,
                         bool revBeam, bool revTrack):
    LinearMapper(beamline, reference, revBeam, revTrack)
{}


IdealMapper::~IdealMapper()
{}


void IdealMapper::getMatrix(FMatrix<double, 6, 6> &matrix) const {
    matrix = itsMap.linearTerms();
}


void IdealMapper::setMatrix(const FMatrix<double, 6, 6> &matrix) {
    itsMap = LinearMap<double, 6>(matrix);
}


void IdealMapper::visitCorrector(const Corrector &corr) {
    // Ignore any kick for ideal orbit.
    applyDrift(flip_s * corr.getElementLength());
}


void IdealMapper::visitPatch(const Patch &/*patch*/) {
    // Ignore patch for ideal orbit.
}


void IdealMapper::visitSeparator(const Separator &sep) {
    // Ignore any kicks for ideal orbit.
    applyDrift(flip_s * sep.getElementLength());
}


void IdealMapper::visitAlignWrapper(const AlignWrapper &wrap) {
    // Ignore misalignments.
    wrap.getElement()->accept(*this);
}


void IdealMapper::visitCorrectorWrapper(const CorrectorWrapper &wrap) {
    // Ignore field errors.
    visitCorrector(wrap.getDesign());
}


void IdealMapper::visitMultipoleWrapper(const MultipoleWrapper &wrap) {
    // Ignore field errors.
    visitMultipole(wrap.getDesign());
}


void IdealMapper::visitRBendWrapper(const RBendWrapper &wrap) {
    // Ignore field errors.
    visitRBend(wrap.getDesign());
}


void IdealMapper::visitSBendWrapper(const SBendWrapper &wrap) {
    // Ignore field errors.
    visitSBend(wrap.getDesign());
}


void IdealMapper::makeFocus(double k, double L, double &c, double &s) {
    double t = k * L * L;
    if(std::abs(t) < 1.0e-4) {
        c = 1.0 - t / 2.0;
        s = L * (1.0 - t / 6.0);
    } else if(k > 0.0) {
        double r = sqrt(k);
        c = cos(r * L);
        s = sin(r * L) / r;
    } else {
        double r = sqrt(- k);
        c = cosh(r * L);
        s = sinh(r * L) / r;
    }
}


void IdealMapper::applyLinearMap
(double length, double kx, double ks, double ky) {
    // Extract the phase coordinates.
    Linear x  = itsMap[X];
    Linear px = itsMap[PX];
    Linear y  = itsMap[Y];
    Linear py = itsMap[PY];

    // Test for zero skew quadrupole component.
    if(ks == 0.0) {
        // Transport coefficients.
        double cx, sx, cy, sy;
        makeFocus(kx, length, cx, sx);
        makeFocus(ky, length, cy, sy);
        double wx = - kx * sx;
        double wy = - ky * sy;

        // Advance through field.
        itsMap[X]  = cx * x + sx * px;
        itsMap[PX] = wx * x + cx * px;
        itsMap[Y]  = cy * y + sy * py;
        itsMap[PY] = wy * y + cy * py;
    } else {
        // Find transformation to principal axes.
        double s1 = (kx + ky) / 2.0;
        double d1 = (kx - ky) / 2.0;
        double root = sqrt(d1 * d1 + ks * ks);
        double c2 = d1 / root;
        double s2 = ks / root;

        // Transport coefficients.
        double cu, su, cv, sv;
        double ku = s1 + (d1 * d1 - ks * ks) / root;
        double kv = s1 - (d1 * d1 - ks * ks) / root;
        makeFocus(ku, length, cu, su);
        makeFocus(kv, length, cv, sv);
        double wu = - ku * su;
        double wv = - kv * sv;

        // Rotate the coordinates to orientation of quadrupole.
        Linear u  = c2 * x  - s2 * y;
        Linear v  = c2 * y  + s2 * x;
        Linear pu = c2 * px - s2 * py;
        Linear pv = c2 * py + s2 * px;

        // Advance through field.
        itsMap[X]  = ((cu + cv) * x + (cu - cv) * u + (su + sv) * px + (su - sv) * pu) / 2.0;
        itsMap[PX] = ((wu + wv) * x + (wu - wv) * u + (cu + cv) * px + (cu - cv) * pu) / 2.0;
        itsMap[Y]  = ((cu + cv) * y - (cu - cv) * v + (su + sv) * py - (su - sv) * pv) / 2.0;
        itsMap[PY] = ((wu + wv) * y - (wu - wv) * v + (cu + cv) * py - (cu - cv) * pv) / 2.0;
    }
}


void IdealMapper::applyMultipoleBody
(double length, double /*refLength*/, const BMultipoleField &field, double /*scale*/) {
    double kx = field.normal(2);
    double ky = - field.normal(2);
    double ks = field.skew(2);
    applyLinearMap(length, kx, ks, ky);
}


void IdealMapper::applySBendBody
(double length, double /*refLength*/, double h,
 const BMultipoleField &field, double /*scale*/) {
    double kx = (h * field.normal(1) + field.normal(2));
    double ky = - field.normal(2);
    double ks = field.skew(2);
    applyLinearMap(length, kx, ks, ky);
}


void IdealMapper::applyThinMultipole
(const BMultipoleField &field, double scale) {
    if(field.order() >= 2) {
        double kn = scale * field.normal(2);
        double ks = scale * field.skew(2);
        itsMap[PX] -= kn * itsMap[X] + ks * itsMap[Y];
        itsMap[PY] -= ks * itsMap[X] - kn * itsMap[Y];
    }
}


void IdealMapper::applyThinSBend
(const BMultipoleField &field, double /*scale*/, double h) {
    double kx = (h * field.normal(1) + field.normal(2));
    double ky = - field.normal(2);
    double ks = field.skew(2);
    itsMap[PX] -= kx * itsMap[X] + ks * itsMap[Y];
    itsMap[PY] -= ks * itsMap[X] + ky * itsMap[Y];
}


void IdealMapper::applyTransform(const Euclid3D &euclid, double /*refLength*/) {
    applyDrift(- euclid.getZ());
    itsMap[PX] += euclid.M(2, 0);
    itsMap[PY] += euclid.M(2, 1);
}
