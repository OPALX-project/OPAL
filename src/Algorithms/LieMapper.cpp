// ------------------------------------------------------------------------
// $RCSfile: LieMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LieMapper
//   The visitor class for building a map of given order for a beamline
//   using a finite-length lenses for all elements.
//   Multipole-like elements are done by Lie series.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 08:16:29 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Algorithms/LieMapper.h"

#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/ElementBase.h"
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

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/DragtFinnMap.h"
#include "FixedAlgebra/FLieGenerator.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "Physics/Physics.h"

class Beamline;
class PartData;
using Physics::c;

typedef FTps<double, 6> Series;


// Class LieMapper
// ------------------------------------------------------------------------

LieMapper::LieMapper(const Beamline &beamline, const PartData &reference,
                     bool backBeam, bool backTrack, int order):
    AbstractMapper(beamline, reference, backBeam, backTrack),
    itsOrder(std::min(order, 6))
{}


LieMapper::~LieMapper()
{}


void LieMapper::getMap(LinearMap<double, 6> &map) const {
    map = itsMap.operator LinearMap<double, 6>();
}


void LieMapper::getMap(FVps<double, 6> &map) const {
    map = itsMap.operator FVps<double, 6>();
}


void LieMapper::getMap(DragtFinnMap<3> &map) const {
    map = itsMap;
}


void LieMapper::setMap(const LinearMap<double, 6> &map) {
    itsMap = DragtFinnMap<3>(map);
}


void LieMapper::setMap(const FVps<double, 6> &map) {
    itsMap = DragtFinnMap<3>(map);
}


void LieMapper::setMap(const DragtFinnMap<3> &map) {
    itsMap = map;
}


void LieMapper::visitBeamBeam(const BeamBeam &map) {
    // *** MISSING *** Map for beam-beam.
}

void LieMapper::visitBeamStripping(const BeamStripping &map) {
    // *** MISSING *** Map for beam stripping.
}

void LieMapper::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}


void LieMapper::visitComponent(const Component &comp) {
    // *** MISSING *** Map for arbitrary component.
}


void LieMapper::visitCorrector(const Corrector &corr) {
    // Get data.
    double length = flip_s * corr.getElementLength();
    double kin = itsReference.getM() / itsReference.getP();
    double scale = (length * itsReference.getQ()) / itsReference.getP();
    BDipoleField field = corr.getField();

    // First order terms (kicks).
    FLieGenerator<double, 3> f_1(1);
    f_1[1] = - scale * field.getBy();
    f_1[2] = - f_1[1] * (length / 2.0);
    f_1[3] =   scale * field.getBx();
    f_1[4] = - f_1[3] * (length / 2.0);

    // Apply map; higher terms ignored.
    DragtFinnMap<3> theMap;
    theMap.assign(f_1);
    theMap.getMatrix()(0, 1) = length;
    theMap.getMatrix()(2, 3) = length;
    theMap.getMatrix()(4, 5) = length * kin * kin;
    itsMap = itsMap.catenate(theMap);
}


void LieMapper::visitDegrader(const Degrader &deg) {
    // The diagnostic has no effect on the map.
    applyDrift(flip_s * deg.getElementLength());
}

void LieMapper::visitDiagnostic(const Diagnostic &diag) {
    // The diagnostic has no effect on the map.
    applyDrift(flip_s * diag.getElementLength());
}


void LieMapper::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void LieMapper::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void LieMapper::visitLambertson(const Lambertson &lamb) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void LieMapper::visitMarker(const Marker &marker) {
    // Do nothing.
}


void LieMapper::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void LieMapper::visitMultipole(const Multipole &mult) {
    double length = mult.getElementLength() * flip_s;
    const BMultipoleField &field = mult.getField();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();

    if(length) {
        // Normal case: Finite-length multipole, field coefficients are B.
        // Apply entrance fringe field.
        applyEntranceFringe(0.0, 0.0, field, scale);

        // Kinematic terms.
        double kin = itsReference.getM() / itsReference.getP();

        // Vector potential in straight reference.
        Series As = buildMultipoleVectorPotential(field) * scale;

        // Build Hamiltonian and factored map, catenate with previous map.
        static const Series px  = Series::makeVariable(AbstractMapper::PX);
        static const Series py  = Series::makeVariable(AbstractMapper::PY);
        static const Series pt1 = Series::makeVariable(AbstractMapper::PT) + 1.0;
        Series pz  = sqrt(pt1 * pt1 - px * px - py * py);
        Series E = sqrt(pt1 * pt1 + kin * kin) / itsReference.getBeta();
        Series H = length * (As + E - pz);
        DragtFinnMap<3> theMap(H);
        itsMap = itsMap.catenate(theMap);

        // Apply exit fringe field.
        applyExitFringe(0.0, 0.0, field, scale);
    } else {
        // Special case: Thin multipole, field coefficients are integral(B*dl).
        // Hamiltonian for thin Multipole is just A_s.
        Series H = buildMultipoleVectorPotential(field) * scale;
        DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H);
        itsMap = itsMap.catenate(theMap);
    }
}


void LieMapper::visitPatch(const Patch &patch) {
    applyTransform(patch.getPatch());
}

void LieMapper::visitProbe(const Probe &prob) {
    // Do nothing.
}


void LieMapper::visitRBend(const RBend &bend) {
    // Geometry.
    const RBendGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();

    if(length) {
        // Finite-length RBend.
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

        // Kinematic terms.
        double kin = itsReference.getM() / itsReference.getP();

        // Vector potential in straight reference.
        Series As = buildMultipoleVectorPotential(field) * scale;

        // Finish Hamiltonian, build factored map, and catenate with previous map.
        static const Series px  = Series::makeVariable(AbstractMapper::PX);
        static const Series py  = Series::makeVariable(AbstractMapper::PY);
        static const Series pt1 = Series::makeVariable(AbstractMapper::PT) + 1.0;
        Series pz  = sqrt(pt1 * pt1 - px * px - py * py);
        Series E  = sqrt(pt1 * pt1 + kin * kin) / itsReference.getBeta();
        Series H = length * (As + E - pz);
        DragtFinnMap<3> theMap(H);
        itsMap = itsMap.catenate(theMap);

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
    } else {
        // Thin RBend.
        double half_angle = flip_s * geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);

        // Transform from in-plane to mid-plane.
        applyTransform(rotat, 0.0);

        // Hamiltonian for thin RBend is just vector potential.
        Series H = buildMultipoleVectorPotential(field) * scale;
        DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H);
        itsMap = itsMap.catenate(theMap);

        // Transform from mid-plane to out-plane.
        applyTransform(rotat, 0.0);
    }
}


void LieMapper::visitParallelPlate(const ParallelPlate &pplate) {
    // Do nothing.
}

void LieMapper::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing here.
}

void LieMapper::visitRFCavity(const RFCavity &as) {
    // Drift through half length.
    double length = flip_s * as.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Get parameters.
    double freq = as.getFrequency();
    double peak = flip_s * as.getAmplitude() / itsReference.getP();

    // Compute Hamiltonian.
    static const Series t = Series::makeVariable(AbstractMapper::T);
    Series H = peak * cos(as.getPhase() + (freq / c) * t);

    // Build map.
    DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H);
    itsMap = itsMap.catenate(theMap);

    // Drift through half length.
    if(length) applyDrift(length / 2.0);
}


void LieMapper::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Map for RF quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void LieMapper::visitSBend(const SBend &bend) {
    const PlanarArcGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();

    if(length) {
        // Apply entrance fringe field.
        applyEntranceFringe(bend.getEntryFaceRotation(),
                            bend.getEntryFaceCurvature(), field, scale);

        // Build Hamiltonian in local coordinates; substitution done later.
        double h = geometry.getCurvature();

        // Kinematic terms.
        double kin = itsReference.getM() / itsReference.getP();

        // Vector potential  (1 + h*x) * As  in curved reference.
        Series As = buildSBendVectorPotential(field, h) * scale;

        // Finish Hamiltonian, build factored map,
        //         and catenate with the previous map.
        static const Series x   = Series::makeVariable(AbstractMapper::X);
        static const Series px  = Series::makeVariable(AbstractMapper::PX);
        static const Series py  = Series::makeVariable(AbstractMapper::PY);
        static const Series pt1 = Series::makeVariable(AbstractMapper::PT) + 1.0;
        Series pz = sqrt(pt1 * pt1 - px * px - py * py);
        Series E  = sqrt(pt1 * pt1 + kin * kin) / itsReference.getBeta();
        Series H  = length * (As + E - (1.0 + h * x) * pz);
        DragtFinnMap<3> theMap(H);
        itsMap = itsMap.catenate(theMap);

        // Apply exit fringe field.
        applyExitFringe(bend.getExitFaceRotation(),
                        bend.getExitFaceCurvature(), field, scale);
    } else {
        double half_angle = geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);

        // Transform from in-plane to mid-plane.
        applyTransform(rotat, 0.0);

        // Hamiltonian for thin SBend is just vector potential,
        // but curvature is non known.
        Series H = buildSBendVectorPotential(field, 0.0) * scale;
        DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H);
        itsMap = itsMap.catenate(theMap);

        // Transform from mid-plane to out-plane.
        applyTransform(rotat, 0.0);
    }
}


void LieMapper::visitSeparator(const Separator &sep) {
    // Get data.
    double length = flip_s * sep.getElementLength();
    double kin = itsReference.getM() / itsReference.getP();
    double scale = (length * itsReference.getQ()) / itsReference.getP();

    // First-order terms (kicks).
    FLieGenerator<double, 3> f_1(1);
    f_1[1] = scale * sep.getEx();
    f_1[2] = - f_1[1] * (length / 2.0);
    f_1[3] = scale * sep.getEy();
    f_1[4] = - f_1[3] * (length / 2.0);

    // Apply map; higher-order terms ignored.
    DragtFinnMap<3> theMap;
    theMap.assign(f_1);
    theMap.getMatrix()(0, 1) = length;
    theMap.getMatrix()(2, 3) = length;
    theMap.getMatrix()(4, 5) = length * kin * kin;
    itsMap = itsMap.catenate(theMap);
}


void LieMapper::visitSeptum(const Septum &sept) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void LieMapper::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double kin = itsReference.getM() / itsReference.getP();
            static const Series px  = Series::makeVariable(AbstractMapper::PX);
            static const Series py  = Series::makeVariable(AbstractMapper::PY);
            static const Series pt1 = Series::makeVariable(AbstractMapper::PT) + 1.0;
            Series pz = sqrt(pt1 * pt1 - px * px - py * py);
            Series H  = length * (sqrt(pt1 * pt1 + kin * kin) - pz);
            DragtFinnMap<3> theMap = DragtFinnMap<3>::factorSimple(H);
            itsMap = itsMap.catenate(theMap);
        } else {
            applyDrift(length);
        }
    }
}


void LieMapper::applyDrift(double length)

{
    // Fetch data.
    double kin = itsReference.getM() / itsReference.getP();

    // Build matrix.
    DragtFinnMap<3> theMap;
    static const Series px  = Series::makeVariable(AbstractMapper::PX);
    static const Series py  = Series::makeVariable(AbstractMapper::PY);
    static const Series pt1 = Series::makeVariable(AbstractMapper::PT) + 1.0;
    Series pz = sqrt(pt1 * pt1 - px * px - py * py);
    Series H  = length * (sqrt(pt1 * pt1 + kin * kin) - pz);
    theMap = DragtFinnMap<3>::factorSimple(H.filter(2, itsOrder));

    // Catenate with previous map.
    itsMap = itsMap.catenate(theMap);
}


void LieMapper::applyEntranceFringe(double angle, double curve,
                                    const BMultipoleField &field,
                                    double scale) {
    // *** MISSING *** map terms for entrance fringe.
    //  double hx = scale * field.normal(1);
    //  double ex = hx * tan(angle);
    //  double ey = hx * tan(angle + itsMap[PX][0]);
}


void LieMapper::applyExitFringe(double angle, double curve,
                                const BMultipoleField &field,
                                double scale) {
    // *** MISSING *** map terms for exit fringe.
    //  double hx = scale * field.normal(1);
    //  double ex = hx * tan(angle);
    //  double ey = hx * tan(angle - itsMap[PX][0]);
}


void LieMapper::applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        // build rotation matrix and compute additional drift length.
        double s2 = (euclid.M(0, 2) * euclid.getX() +
                     euclid.M(1, 2) * euclid.getY() +
                     euclid.M(2, 2) * euclid.getZ()) / euclid.M(2, 2);
        double kin = itsReference.getM() / itsReference.getP();

        // f1 terms (kicks).
        FLieGenerator<double, 3> f_1(1);
        f_1[1] = euclid.M(0, 2) / euclid.M(2, 2);
        f_1[2] = - euclid.getX() + f_1[1] * s2;
        f_1[3] = euclid.M(1, 2) / euclid.M(2, 2);
        f_1[4] = - euclid.getY() + f_1[3] * s2;
        f_1[5] = 0.0;
        f_1[6] = euclid.getZ();

        // Transfer matrix.
        FMatrix<double, 6, 6> F;
        F(1, 1) = euclid.M(0, 0);
        F(1, 3) = euclid.M(1, 0);
        F(1, 5) = euclid.M(2, 0);
        F(3, 1) = euclid.M(0, 1);
        F(3, 3) = euclid.M(1, 1);
        F(3, 5) = euclid.M(2, 1);
        F(5, 5) = 1.0;

        F(0, 0) = euclid.M(1, 1) / euclid.M(2, 2);
        F(0, 2) = F(0, 0) * s2;
        F(0, 2) = - euclid.M(0, 1) / euclid.M(2, 2);
        F(0, 3) = F(0, 2) * s2;
        F(2, 0) = - euclid.M(1, 0) / euclid.M(2, 2);
        F(2, 1) = F(2, 0) * s2;
        F(2, 2) = euclid.M(0, 0) / euclid.M(2, 2);
        F(2, 3) = F(2, 2) * s2;
        F(4, 0) = euclid.M(0, 2) / euclid.M(2, 2);
        F(4, 1) = F(4, 0) * s2;
        F(4, 2) = euclid.M(1, 2) / euclid.M(2, 2);
        F(4, 3) = F(4, 2) * s2;
        F(4, 4) = 1.0;
        F(4, 5) = - s2 * kin * kin;

        // Store generator and matrix.
        DragtFinnMap <3> theMap;
        theMap.assign(f_1);
        theMap.assign(F);
        itsMap = itsMap.catenate(theMap);
    }
}