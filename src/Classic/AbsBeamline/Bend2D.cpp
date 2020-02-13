// ------------------------------------------------------------------------
// $RCSfile: Bend2D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: Bend2D
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Bend2D.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include "Structure/MeshGenerator.h"

#include "gsl/gsl_poly.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class Bend2D
// ------------------------------------------------------------------------

Bend2D::Bend2D():
    Bend2D("")
{}

Bend2D::Bend2D(const Bend2D &right):
    BendBase(right),
    messageHeader_m(" * "),
    pusher_m(right.pusher_m),
    fieldmap_m(right.fieldmap_m),
    designRadius_m(right.designRadius_m),
    exitAngle_m(right.exitAngle_m),
    fieldIndex_m(right.fieldIndex_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    widthEntranceFringe_m(right.widthEntranceFringe_m),
    widthExitFringe_m(right.widthExitFringe_m),
    reinitialize_m(right.reinitialize_m),
    entranceParameter1_m(right.entranceParameter1_m),
    entranceParameter2_m(right.entranceParameter2_m),
    entranceParameter3_m(right.entranceParameter3_m),
    exitParameter1_m(right.exitParameter1_m),
    exitParameter2_m(right.exitParameter2_m),
    exitParameter3_m(right.exitParameter3_m),
    engeCoeffsEntry_m(right.engeCoeffsEntry_m),
    engeCoeffsExit_m(right.engeCoeffsExit_m),
    entryFieldValues_m(NULL),
    exitFieldValues_m(NULL),
    entryFieldAccel_m(NULL),
    exitFieldAccel_m(NULL),
    deltaBeginEntry_m(right.deltaBeginEntry_m),
    deltaEndEntry_m(right.deltaEndEntry_m),
    polyOrderEntry_m(right.polyOrderEntry_m),
    deltaBeginExit_m(right.deltaBeginExit_m),
    deltaEndExit_m(right.deltaEndExit_m),
    polyOrderExit_m(right.polyOrderExit_m),
    cosEntranceAngle_m(right.cosEntranceAngle_m),
    sinEntranceAngle_m(right.sinEntranceAngle_m),
    tanEntranceAngle_m(right.tanEntranceAngle_m),
    tanExitAngle_m(right.tanExitAngle_m),
    beginToEnd_m(right.beginToEnd_m),
    beginToEnd_lcs_m(right.beginToEnd_lcs_m),
    toEntranceRegion_m(right.toEntranceRegion_m),
    toExitRegion_m(right.toExitRegion_m),
    computeAngleTrafo_m(right.computeAngleTrafo_m),
    maxAngle_m(right.maxAngle_m),
    nSlices_m(right.nSlices_m){

    setElType(isDipole);

}

Bend2D::Bend2D(const std::string &name):
    BendBase(name),
    messageHeader_m(" * "),
    pusher_m(),
    fieldmap_m(NULL),
    designRadius_m(0.0),
    exitAngle_m(0.0),
    fieldIndex_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    widthEntranceFringe_m(0.0),
    widthExitFringe_m(0.0),
    reinitialize_m(false),
    entranceParameter1_m(0.0),
    entranceParameter2_m(0.0),
    entranceParameter3_m(0.0),
    exitParameter1_m(0.0),
    exitParameter2_m(0.0),
    exitParameter3_m(0.0),
    entryFieldValues_m(NULL),
    exitFieldValues_m(NULL),
    entryFieldAccel_m(NULL),
    exitFieldAccel_m(NULL),
    deltaBeginEntry_m(0.0),
    deltaEndEntry_m(0.0),
    polyOrderEntry_m(0),
    deltaBeginExit_m(0.0),
    deltaEndExit_m(0.0),
    polyOrderExit_m(0),
    cosEntranceAngle_m(1.0),
    sinEntranceAngle_m(0.0),
    tanEntranceAngle_m(0.0),
    tanExitAngle_m(0.0),
    maxAngle_m(0.0),
    nSlices_m(1){

    setElType(isDipole);

}

Bend2D::~Bend2D() {
    if (entryFieldAccel_m != NULL) {
        for (unsigned int i = 0; i < 3u; ++ i) {
            gsl_spline_free(entryFieldValues_m[i]);
            gsl_spline_free(exitFieldValues_m[i]);
            entryFieldValues_m[i] = NULL;
            exitFieldValues_m[i] = NULL;
        }
        delete[] entryFieldValues_m;
        delete[] exitFieldValues_m;

        gsl_interp_accel_free(entryFieldAccel_m);
        gsl_interp_accel_free(exitFieldAccel_m);
    }
}
/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */
bool Bend2D::apply(const size_t &i,
                 const double &t,
                 Vector_t &E,
                 Vector_t &B) {
    Vector_t P;
    return apply(RefPartBunch_m->R[i], P, t, E, B);
}

bool Bend2D::apply(const Vector_t &R,
                   const Vector_t &/*P*/,
                   const double &/*t*/,
                   Vector_t &/*E*/,
                   Vector_t &B) {

    if(designRadius_m > 0.0) {

        Vector_t X = R;
        Vector_t bField(0.0);
        if (!calculateMapField(X, bField)) {

            B += fieldAmplitude_m * bField;

            return false;
        }

        return true;
    }
    return false;

}

bool Bend2D::applyToReferenceParticle(const Vector_t &R,
                                    const Vector_t &P,
                                    const double &t,
                                    Vector_t &E,
                                    Vector_t &B) {
    return apply(R, P, t, E, B);
}

void Bend2D::goOnline(const double &) {
    online_m = true;
}

void Bend2D::initialise(PartBunchBase<double, 3> *bunch,
                      double &startField,
                      double &endField) {

    Inform msg(messageHeader_m.c_str(), *gmsg);

    if(initializeFieldMap(msg)) {
        std::string name = Util::toUpper(getName()) + " ";
        msg << level2
            << "======================================================================\n"
            << "===== " << std::left << std::setw(64) << std::setfill('=') << name << "\n"
            << "======================================================================\n";

        setupPusher(bunch);
        readFieldMap(msg);
        setupBendGeometry(msg, startField, endField);

        double bendAngleX = 0.0;
        double bendAngleY = 0.0;
        calculateRefTrajectory(bendAngleX, bendAngleY);
        print(msg, bendAngleX, bendAngleY);

        // Pass start and end of field to calling function.
        startField = startField_m;
        endField = endField_m;

        startField_m -= elementEdge_m;
        endField_m -= elementEdge_m;
        elementEdge_m = 0.0;

        setupFringeWidths();
        msg << level2
            << "======================================================================\n"
            << "======================================================================\n"
            << endl;
    } else {
        ERRORMSG("There is something wrong with your field map \""
                 << fileName_m
                 << "\"");
        endField = startField - 1.0e-3;
    }
}

void Bend2D::adjustFringeFields(double ratio) {
    findChordLength(chordLength_m);

    double delta = std::abs(entranceParameter1_m - entranceParameter2_m);
    entranceParameter1_m = entranceParameter2_m - delta * ratio;

    delta = std::abs(entranceParameter2_m - entranceParameter3_m);
    entranceParameter3_m = entranceParameter2_m + delta * ratio;

    delta = std::abs(exitParameter1_m - exitParameter2_m);
    exitParameter1_m = exitParameter2_m - delta * ratio;

    delta = std::abs(exitParameter2_m - exitParameter3_m);
    exitParameter3_m = exitParameter2_m + delta * ratio;

    setupFringeWidths();
}

double Bend2D::calculateBendAngle() {

    const double mass = RefPartBunch_m->getM();
    const double gamma = designEnergy_m / mass + 1.0;
    const double betaGamma = sqrt(pow(gamma, 2.0) - 1.0);
    // const double beta = betaGamma / gamma;
    const double deltaT = RefPartBunch_m->getdT();
    const double cdt = Physics::c * deltaT;
    // Integrate through field for initial angle.
    Vector_t oldX;
    Vector_t X = -deltaBeginEntry_m * Vector_t(tan(entranceAngle_m), 0.0, 1.0);
    Vector_t P = betaGamma * Vector_t(sin(entranceAngle_m), 0.0, cos(entranceAngle_m));
    double deltaS = 0.0;
    double bendLength = endField_m - startField_m;
    const Vector_t eField(0.0);

    while(deltaS < bendLength) {
        oldX = X;
        X /= cdt;
        pusher_m.push(X, P, deltaT);
        X *= cdt;

        Vector_t bField(0.0, 0.0, 0.0);
        calculateMapField(X, bField);
        bField = fieldAmplitude_m * bField;

        X /= cdt;
        pusher_m.kick(X, P, eField, bField, deltaT);

        pusher_m.push(X, P, deltaT);
        X *= cdt;

        deltaS += euclidean_norm(X - oldX);

    }

    double angle =  -atan2(P(0), P(2)) + entranceAngle_m;

    return angle;

}

void Bend2D::calcEngeFunction(double zNormalized,
                             const std::vector<double> &engeCoeff,
                             int polyOrder,
                             double &engeFunc,
                             double &engeFuncDeriv,
                             double &engeFuncSecDerivNorm) {

    double expSum = 0.0;
    double expSumDeriv = 0.0;
    double expSumSecDeriv = 0.0;

    if(polyOrder >= 2) {

        expSum = engeCoeff.at(0)
            + engeCoeff.at(1) * zNormalized;
        expSumDeriv = engeCoeff.at(1);

        for(int index = 2; index <= polyOrder; index++) {
            expSum += engeCoeff.at(index) * pow(zNormalized, index);
            expSumDeriv += index * engeCoeff.at(index)
                * pow(zNormalized, index - 1);
            expSumSecDeriv += index * (index - 1) * engeCoeff.at(index)
                * pow(zNormalized, index - 2);
        }

    } else if(polyOrder == 1) {

        expSum = engeCoeff.at(0)
            + engeCoeff.at(1) * zNormalized;
        expSumDeriv = engeCoeff.at(1);

    } else
        expSum = engeCoeff.at(0);

    double engeExp = exp(expSum);
    engeFunc = 1.0 / (1.0 + engeExp);

    if(!std::isnan(engeFunc)) {

        expSumDeriv /= gap_m;
        expSumSecDeriv /= pow(gap_m, 2.0);
        double engeExpDeriv = expSumDeriv * engeExp;
        double engeExpSecDeriv = (expSumSecDeriv + pow(expSumDeriv, 2.0))
            * engeExp;
        double engeFuncSq = pow(engeFunc, 2.0);

        engeFuncDeriv = -engeExpDeriv * engeFuncSq;
        if (std::isnan(engeFuncDeriv) || std::isinf(engeFuncDeriv))
            engeFuncDeriv = 0.0;

        engeFuncSecDerivNorm = -engeExpSecDeriv * engeFunc
            + 2.0 * pow(engeExpDeriv, 2.0)
            * engeFuncSq;
        if (std::isnan(engeFuncSecDerivNorm) || std::isinf(engeFuncSecDerivNorm))
            engeFuncSecDerivNorm = 0.0;

    } else {
        engeFunc = 0.0;
        engeFuncDeriv = 0.0;
        engeFuncSecDerivNorm = 0.0;

    }
}

// :FIXME: is this correct?
Vector_t Bend2D::calcCentralField(const Vector_t &/*R*/, double /*deltaX*/) {

    Vector_t B(0, 0, 0);
    //double nOverRho = fieldIndex_m / designRadius_m;
    //double expFactor = exp(-nOverRho * deltaX);
    //double bxBzFactor = expFactor * nOverRho * R(1);
    //Vector_t rotationCenter(-designRadius_m, R(1), 0.0);
    //double cosangle = dot(R - rotationCenter, Vector_t(1, 0, 0)) / euclidean_norm(R - rotationCenter);

    //B(0) = -bxBzFactor * cosangle;
    //B(1) = expFactor * (1.0 - pow(nOverRho * R(1), 2.0) / 2.0);
    //B(2) = -bxBzFactor * sqrt(1 - std::pow(cosangle, 2));

    B(1) = 1.0;

    return B;
}

Vector_t Bend2D::calcEntranceFringeField(const Vector_t &R,
                                         double /*deltaX*/) {

    const CoordinateSystemTrafo toEntranceRegion(Vector_t(0, 0, entranceParameter2_m),
                                                 Quaternion(0, 0, 1, 0));
    const Vector_t Rprime = toEntranceRegion.transformTo(R);

    Vector_t B(0.0);

    if (Rprime(2) <= -entranceParameter3_m) {
        B(1) = 1;
    } else if (Rprime(2) < -entranceParameter1_m) {
        double engeFunc = gsl_spline_eval(entryFieldValues_m[0], Rprime(2), entryFieldAccel_m);
        double engeFuncDeriv = gsl_spline_eval(entryFieldValues_m[1], Rprime(2), entryFieldAccel_m);
        double engeFuncSecDeriv = gsl_spline_eval(entryFieldValues_m[2], Rprime(2), entryFieldAccel_m);

        //double nOverRho = fieldIndex_m / designRadius_m;
        //double expFactor = exp(-nOverRho * deltaX);
        //double trigFactor = pow(nOverRho, 2.0) + engeFuncSecDerivNorm;

        //double bXEntrance = -nOverRho * expFactor * Rprime(1) * engeFunc;
        //double bYEntrance = (expFactor * engeFunc *
        //                      (1.0  - 0.5 * trigFactor * pow(Rprime(1), 2.0)));
        //double bZEntrance = expFactor * Rprime(1) * engeFuncDeriv;


        // B(1) = (engeFunc *
        //  (1.0 - 0.5 * engeFuncSecDerivNorm * pow(Rprime(1), 2.0)));

        B(1) = (engeFunc - 0.5 * engeFuncSecDeriv * pow(Rprime(1), 2.0));

        B(2) = engeFuncDeriv * Rprime(1);
    }

    return toEntranceRegion.rotateFrom(B);
}

Vector_t Bend2D::calcExitFringeField(const Vector_t &R,
                                     double /*deltaX*/) {

    const CoordinateSystemTrafo fromEndToExitRegion(Vector_t(0, 0, exitParameter2_m),
                                                    Quaternion(1, 0, 0, 0));
    const CoordinateSystemTrafo toExitRegion = (fromEndToExitRegion *
                                                 getBeginToEnd_local());
    const Vector_t Rprime = toExitRegion.transformTo(R);

    Vector_t B(0.0);

    if (Rprime(2) <= exitParameter1_m) {
        B(1) = 1;
    } else if (Rprime(2) < exitParameter3_m) {
        double engeFunc = gsl_spline_eval(exitFieldValues_m[0], Rprime(2), exitFieldAccel_m);
        double engeFuncDeriv = gsl_spline_eval(exitFieldValues_m[1], Rprime(2), exitFieldAccel_m);
        double engeFuncSecDeriv = gsl_spline_eval(exitFieldValues_m[2], Rprime(2), exitFieldAccel_m);

        //double nOverRho = fieldIndex_m / designRadius_m;
        //double expFactor = exp(-nOverRho * deltaX);
        //double trigFactor = pow(nOverRho, 2.0) + engeFuncSecDerivNorm;

        //double bXExit = -nOverRho * expFactor * Rprime(1) * engeFunc;
        //double bYExit = (expFactor * engeFunc *
        //                 (1.0 - 0.5 * trigFactor * pow(Rprime(1), 2.0)));
        //double bZExit = expFactor * Rprime(1) * engeFuncDeriv;

        //B(1) = (engeFunc *
        //        (1.0 - 0.5 * engeFuncSecDerivNorm * pow(Rprime(1), 2.0)));
        B(1) = (engeFunc - 0.5 * engeFuncSecDeriv * pow(Rprime(1), 2.0));

        B(2) = engeFuncDeriv * Rprime(1);
    }
    return toExitRegion.rotateFrom(B);
}

bool Bend2D::calculateMapField(const Vector_t &R, Vector_t &B) {

    B = Vector_t(0.0);
    bool verticallyInside = (std::abs(R(1)) < 0.5 * gap_m);
    bool horizontallyInside = false;
    Vector_t rotationCenter(-designRadius_m * cosEntranceAngle_m, R(1), designRadius_m * sinEntranceAngle_m);
    if(inMagnetCentralRegion(R)) {
        if (verticallyInside) {
            double deltaX = 0.0;//euclidean_norm(R - rotationCenter) - designRadius_m;
            bool inEntranceRegion = isPositionInEntranceField(R);
            bool inExitRegion = isPositionInExitField(R);

            if (inEntranceRegion && inExitRegion) {
                Vector_t Rp = transformToEntranceRegion(R);
                Vector_t Rpp = transformToExitRegion(R);

                if (std::abs(Rp[2]) < std::abs(Rpp[2])) {
                    inExitRegion = false;
                } else {
                    inEntranceRegion = false;
                }
            }
            if (inEntranceRegion) {
                B = calcEntranceFringeField(R, deltaX);
            } else if (inExitRegion) {
                B = calcExitFringeField(R, deltaX);
            } else {
                B = calcCentralField(R, deltaX);
            }

            return false;
        }
        return true;
    }

    Vector_t BEntrance(0.0), BExit(0.0);
    verticallyInside = (std::abs(R(1)) < gap_m);

    bool inEntranceRegion = inMagnetEntranceRegion(R);
    bool inExitRegion = inMagnetExitRegion(R);

    if (inEntranceRegion) {
        horizontallyInside = true;
        if (verticallyInside) {
            BEntrance = calcEntranceFringeField(R, R(0));
        }
    }

    if (inExitRegion) {
        horizontallyInside = true;
        if (verticallyInside) {
            Vector_t Rprime = getBeginToEnd_local().transformTo(R);

            BExit = calcExitFringeField(R, Rprime(0));
        }
    }

    B = BEntrance + BExit;

    bool hitMaterial = (horizontallyInside && (!verticallyInside));
    return hitMaterial;
}

void Bend2D::calculateRefTrajectory(double &angleX, double &/*angleY*/) {

    const double mass = RefPartBunch_m->getM();
    const double gamma = designEnergy_m / mass + 1.;
    const double betaGamma = sqrt(gamma * gamma - 1.);
    const double dt = RefPartBunch_m->getdT();

    std::ofstream trajectoryOutput;
    if (Options::writeBendTrajectories && Ippl::myNode() == 0) {
        trajectoryOutput.open("data/" + OpalData::getInstance()->getInputBasename() + "_" + getName() + "_traj.dat");
        trajectoryOutput.precision(12);
        trajectoryOutput << "# " << std::setw(18) << "s"
                         << std::setw(20) << "x"
                         << std::setw(20) << "z"
                         << std::setw(20) << "By"
                         << std::endl;
    }

    double zRotation = rotationZAxis_m;
    Quaternion toStandard = Quaternion(cos(0.5 * zRotation), sin(0.5 * zRotation) * Vector_t(0, 0, -1));

    Vector_t X = -deltaBeginEntry_m * Vector_t(tan(entranceAngle_m), 0.0, 1.0);
    Vector_t P = betaGamma * Vector_t(sin(entranceAngle_m), 0.0, cos(entranceAngle_m));

    if(!refTrajMap_m.empty())
        refTrajMap_m.clear();

    refTrajMap_m.push_back(X);

    const Vector_t eField(0.0);
    const double stepSize = betaGamma / gamma * Physics::c * dt;
    const double bendLength = endField_m - startField_m;
    double deltaS = 0.0;
    while(deltaS < bendLength) {

        X /= Vector_t(Physics::c * dt);
        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        Vector_t bField(0.0, 0.0, 0.0);
        Vector_t XInBendFrame = X;

        calculateMapField(XInBendFrame, bField);
        bField = fieldAmplitude_m * bField;

        if (Options::writeBendTrajectories && Ippl::myNode() == 0) {
            trajectoryOutput << std::setw(20) << deltaS + 0.5 * stepSize
                             << std::setw(20) << X(0)
                             << std::setw(20) << X(2)
                             << std::setw(20) << bField(1)
                             << std::endl;
        }

        X /= Vector_t(Physics::c * dt);
        pusher_m.kick(X, P, eField, bField, dt);

        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        refTrajMap_m.push_back(X);

        deltaS += stepSize;
    }

    angleX = -atan2(P(0), P(2)) + entranceAngle_m;
}

double Bend2D::estimateFieldAdjustmentStep(double actualBendAngle,
                                          double mass,
                                          double betaGamma) {

    // Estimate field adjustment step.
    double effectiveLength = angle_m * designRadius_m;
    double tmp1 = betaGamma * mass / (2.0 * effectiveLength * Physics::c);
    double tmp2 = pow(fieldAmplitude_m / tmp1, 2.0);
    double fieldStep = (angle_m - actualBendAngle) * tmp1;
    if (tmp2 < 1.0) {
        fieldStep = (angle_m - actualBendAngle) * tmp1 * std::sqrt(1.0 - tmp2);
    }
    fieldStep *= fieldAmplitude_m / std::abs(fieldAmplitude_m);

    return fieldStep;

}

void Bend2D::findBendEffectiveLength(double startField, double endField) {

    /*
     * Use an iterative procedure to set the width of the
     * default field map for the defined field amplitude
     * and bend angle.
     */
    setEngeOriginDelta(0.0);
    setFieldCalcParam();
    setFieldBoundaries(startField, endField);

    double actualBendAngle = calculateBendAngle();
    double error = std::abs(actualBendAngle - angle_m);

    if(error > 1.0e-6) {

        double deltaStep = 0.0;
        if(std::abs(actualBendAngle) < std::abs(angle_m))
            deltaStep = -gap_m / 2.0;
        else
            deltaStep = gap_m / 2.0;

        double delta1 = 0.0;
        double bendAngle1 = actualBendAngle;

        double delta2 = deltaStep;
        setEngeOriginDelta(delta2);
        setFieldCalcParam();
        setFieldBoundaries(startField, endField);
        double bendAngle2 = calculateBendAngle();

        if(std::abs(bendAngle1) > std::abs(angle_m)) {
            while(std::abs(bendAngle2) > std::abs(angle_m)) {
                delta2 += deltaStep;
                setEngeOriginDelta(delta2);
                setFieldCalcParam();
                setFieldBoundaries(startField, endField);
                bendAngle2 = calculateBendAngle();
            }
        } else {
            while(std::abs(bendAngle2) < std::abs(angle_m)) {
                delta2 += deltaStep;
                setEngeOriginDelta(delta2);
                setFieldCalcParam();
                setFieldBoundaries(startField, endField);
                bendAngle2 = calculateBendAngle();
            }
        }

        // Now we should have the proper field map width bracketed.
        unsigned int iterations = 1;
        error = std::abs(actualBendAngle - angle_m);
        while(error > 1.0e-6 && iterations < 100) {

            double delta = (delta1 + delta2) / 2.0;
            setEngeOriginDelta(delta);
            setFieldCalcParam();
            setFieldBoundaries(startField, endField);
            double newBendAngle = calculateBendAngle();

            error = std::abs(newBendAngle - angle_m);

            if(error > 1.0e-6) {

                if(bendAngle1 - angle_m < 0.0) {

                    if(newBendAngle - angle_m < 0.0) {
                        bendAngle1 = newBendAngle;
                        delta1 = delta;
                    } else {
                        // bendAngle2 = newBendAngle;
                        delta2 = delta;
                    }

                } else {

                    if(newBendAngle - angle_m < 0.0) {
                        // bendAngle2 = newBendAngle;
                        delta2 = delta;
                    } else {
                        bendAngle1 = newBendAngle;
                        delta1 = delta;
                    }
                }
            }
            iterations++;
        }
    }
}

void Bend2D::findBendStrength(double mass,
                              double /*gamma*/,
                              double betaGamma,
                              double /*charge*/) {

    /*
     * Use an iterative procedure to set the magnet field amplitude
     * for the defined bend angle.
     */
    const double tolerance = 1e-7;
    double actualBendAngle = calculateBendAngle();
    double error = std::abs(actualBendAngle - angle_m) * Physics::rad2deg;
    if (error < tolerance)
        return;

    double fieldStep = estimateFieldAdjustmentStep(actualBendAngle,
                                                   mass,
                                                   betaGamma);
    double amplitude1 = fieldAmplitude_m;
    double bendAngle1 = actualBendAngle;

    double amplitude2 = fieldAmplitude_m + fieldStep;
    fieldAmplitude_m = amplitude2;
    double bendAngle2 = calculateBendAngle();

    if(std::abs(bendAngle1) > std::abs(angle_m)) {
        while(std::abs(bendAngle2) > std::abs(angle_m)) {
            amplitude1 = amplitude2;
            bendAngle1 = bendAngle2;

            amplitude2 += fieldStep;
            fieldAmplitude_m = amplitude2;
            bendAngle2 = calculateBendAngle();
        }
    } else {
        while(std::abs(bendAngle2) < std::abs(angle_m)) {
            amplitude1 = amplitude2;
            bendAngle1 = bendAngle2;

            amplitude2 += fieldStep;
            fieldAmplitude_m = amplitude2;
            bendAngle2 = calculateBendAngle();
        }
    }

    // Now we should have the proper field amplitude bracketed.
    unsigned int iterations = 1;
    while(error > tolerance && iterations < 100) {

        fieldAmplitude_m = (amplitude1 + amplitude2) / 2.0;
        double newBendAngle = calculateBendAngle();

        error = std::abs(newBendAngle - angle_m) * Physics::rad2deg;

        if(error > tolerance) {

            if(bendAngle1 - angle_m < 0.0) {

                if(newBendAngle - angle_m < 0.0) {
                    bendAngle1 = newBendAngle;
                    amplitude1 = fieldAmplitude_m;
                } else {
                    // bendAngle2 = newBendAngle;
                    amplitude2 = fieldAmplitude_m;
                }

            } else {

                if(newBendAngle - angle_m < 0.0) {
                    // bendAngle2 = newBendAngle;
                    amplitude2 = fieldAmplitude_m;
                } else {
                    bendAngle1 = newBendAngle;
                    amplitude1 = fieldAmplitude_m;
                }
            }
        }
        iterations++;
    }
}

bool Bend2D::findIdealBendParameters(double chordLength) {

    double refMass = RefPartBunch_m->getM();
    double refGamma = designEnergy_m / refMass + 1.0;
    double refBetaGamma = sqrt(pow(refGamma, 2.0) - 1.0);
    double refCharge = RefPartBunch_m->getQ();
    bool reinitialize = false;

    if(angle_m != 0.0) {

        if(angle_m < 0.0) {
            // Negative angle is a positive bend rotated 180 degrees.
            entranceAngle_m = copysign(1, angle_m) * entranceAngle_m;
            setExitAngle(copysign(1, angle_m) * exitAngle_m);
            angle_m = std::abs(angle_m);
            rotationZAxis_m += Physics::pi;
        }
        designRadius_m = chordLength / (2.0 * std::sin(angle_m / 2.0));

        fieldAmplitude_m = ((refCharge / std::abs(refCharge)) *
                            refBetaGamma * refMass /
                            (Physics::c * designRadius_m));
        reinitialize = true;
    } else {

        rotationZAxis_m += atan2(bX_m, bY_m);
        if(refCharge < 0.0) {
            rotationZAxis_m -= Physics::pi;
        }

        fieldAmplitude_m = (refCharge *
                            std::abs(sqrt(pow(bY_m, 2.0) + pow(bX_m, 2.0)) / refCharge));
        designRadius_m = std::abs(refBetaGamma * refMass / (Physics::c * fieldAmplitude_m));
        double bendAngle = 2.0 * std::asin(chordLength / (2.0 * designRadius_m));

        angle_m = bendAngle;

        reinitialize = false;
    }
    return reinitialize;
}

bool Bend2D::initializeFieldMap(Inform &msg) {

    fieldmap_m = Fieldmap::getFieldmap(fileName_m, fast_m);

    if(fieldmap_m != NULL) {
        if(fileName_m != "1DPROFILE1-DEFAULT")
            return true;
        else
            return setupDefaultFieldMap(msg);

    } else
        return false;

}

bool Bend2D::inMagnetCentralRegion(const Vector_t &R) const {

    Vector_t rotationCenter(-designRadius_m * cosEntranceAngle_m, R(1), designRadius_m * sinEntranceAngle_m);
    double distFromRotCenter = euclidean_norm(R - rotationCenter);
    Vector_t Rprime = getBeginToEnd_local().transformTo(R);
    Vector_t Rpprime = computeAngleTrafo_m.transformTo(R);

    double effectiveAngle = fmod(Physics::two_pi - atan2(Rpprime(0), Rpprime(2)), Physics::two_pi);

    if (std::abs(distFromRotCenter - designRadius_m) < 0.5 * aperture_m.second[0] &&
        effectiveAngle >= 0.0 && effectiveAngle < maxAngle_m) {
        if (effectiveAngle < 0.5 * maxAngle_m) return R(2) >= 0.0;
        return Rprime(2) < 0.0;
    }

    return false;
}

bool Bend2D::inMagnetEntranceRegion(const Vector_t &R) const {

    return (R(2) >= entranceParameter1_m &&
            R(2) < 0.0 &&
            std::abs(R(0)) < aperture_m.second[0]);
}

bool Bend2D::inMagnetExitRegion(const Vector_t &R) const {

    Vector_t Rprime = getBeginToEnd_local().transformTo(R);

    return (Rprime(2) >= 0 &&
            Rprime(2) < exitParameter3_m &&
            std::abs(Rprime(0)) < aperture_m.second[0]);
}

bool Bend2D::isPositionInEntranceField(const Vector_t &R) const {

    return (polyOrderEntry_m >= 0 &&
            R(2) >= entranceParameter1_m &&
            R(2) < entranceParameter3_m);
}

bool Bend2D::isPositionInExitField(const Vector_t &R) const {
    Vector_t Rprime = getBeginToEnd_local().transformTo(R);

    return (polyOrderExit_m >= 0 &&
            Rprime(2) >= exitParameter1_m &&
            Rprime(2) < exitParameter3_m);
}

void Bend2D::print(Inform &msg, double bendAngleX, double bendAngleY) {
    msg << level2 << "\n"
        << "Start of field map:      "
        << startField_m
        << " m (in s coordinates)"
        << "\n";
    msg << "End of field map:        "
        << endField_m
        << " m (in s coordinates)"
        << "\n";
    msg << "Entrance edge of magnet: "
        << elementEdge_m
        << " m (in s coordinates)"
        << "\n";
    msg << "\n"
        << "Reference Trajectory Properties"
        << "\n"
        << "======================================================================"
        << "\n\n";
    msg << "Bend angle magnitude:    "
        << angle_m
        << " rad ("
        << angle_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Entrance edge angle:     "
        << entranceAngle_m
        << " rad ("
        << entranceAngle_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Exit edge angle:         "
        << exitAngle_m
        << " rad ("
        << exitAngle_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Bend design radius:      "
        << designRadius_m
        << " m"
        << "\n";
    msg << "Bend design energy:      "
        << designEnergy_m
        << " eV"
        << "\n";
    msg << "\n"
        << "Bend Field and Rotation Properties"
        << "\n"
        << "======================================================================"
        << "\n\n";
    msg << "Field amplitude:         "
        << fieldAmplitude_m
        << " T"
        << "\n";
    msg << "Field index:  "
        << fieldIndex_m
        << "\n";
    msg << "Rotation about z axis:   "
        << rotationZAxis_m
        << " rad ("
        << rotationZAxis_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "\n"
        << "Reference Trajectory Properties Through Bend Magnet with Fringe Fields"
        << "\n"
        << "======================================================================"
        << "\n\n";
    msg << "Reference particle is bent: "
        << bendAngleX
        << " rad ("
        << bendAngleX * 180.0 / Physics::pi
        << " degrees) in x plane"
        << "\n";
    msg << "Reference particle is bent: "
        << bendAngleY
        << " rad ("
        << bendAngleY * 180.0 / Physics::pi
        << " degrees) in y plane"
        << "\n";

    if(fileName_m == "1DPROFILE1-DEFAULT") {
        msg << "\n"
            << "Effective Field Map\n"
            << "======================================================================\n"
            << "\n"
            << "1DProfile1 " << polyOrderEntry_m << " " << polyOrderExit_m << " " << gap_m * 100 << "\n"
            << entranceParameter1_m * 100 << " " << entranceParameter2_m * 100 << " " << entranceParameter3_m * 100 << " 0\n"
            << exitParameter1_m * 100 << " " << exitParameter2_m * 100 << " " << exitParameter3_m * 100 << " 0\n";
        for (auto coef: engeCoeffsEntry_m) {
            msg << coef << "\n";
        }
        for (auto coef: engeCoeffsExit_m) {
            msg << coef << "\n";
        }
    }
    msg << endl;
}


void Bend2D::readFieldMap(Inform &msg) {
    msg << level2 << getName() << " using file ";
    fieldmap_m->getInfo(&msg);

    Fieldmap::readMap(fileName_m);
    fieldmap_m->get1DProfile1EntranceParam(entranceParameter1_m,
                                           entranceParameter2_m,
                                           entranceParameter3_m);
    fieldmap_m->get1DProfile1ExitParam(exitParameter1_m,
                                       exitParameter2_m,
                                       exitParameter3_m);

    setGapFromFieldMap();
    fieldmap_m->get1DProfile1EngeCoeffs(engeCoeffsEntry_m,
                                        engeCoeffsExit_m);
    polyOrderEntry_m = engeCoeffsEntry_m.size() - 1;
    polyOrderExit_m = engeCoeffsExit_m.size() - 1;

    double stepSize = Physics::c * 1e-12;
    double entryLength = std::abs(entranceParameter3_m - entranceParameter1_m);
    unsigned int numStepsEntry = std::ceil(entryLength / stepSize) + 3;
    double stepSizeEntry = entryLength / (numStepsEntry - 3);
    std::vector<double> zvalsEntry(numStepsEntry);
    std::vector<std::vector<double> > fieldValuesEntry(3);

    double exitLength = std::abs(exitParameter3_m - exitParameter1_m);
    unsigned int numStepsExit = std::ceil(exitLength / stepSize) + 3;
    double stepSizeExit = exitLength / (numStepsExit - 3);
    std::vector<double> zvalsExit(numStepsExit);
    std::vector<std::vector<double> > fieldValuesExit(3);

    entryFieldValues_m = new gsl_spline*[3];
    exitFieldValues_m = new gsl_spline*[3];
    entryFieldAccel_m = gsl_interp_accel_alloc();
    exitFieldAccel_m = gsl_interp_accel_alloc();

    for (unsigned int i = 0; i < 3u; ++ i) {
        fieldValuesEntry[i].resize(numStepsEntry);
        fieldValuesExit[i].resize(numStepsExit);

        entryFieldValues_m[i] = gsl_spline_alloc(gsl_interp_cspline, numStepsEntry);
        exitFieldValues_m[i] = gsl_spline_alloc(gsl_interp_cspline, numStepsExit);
    }

    for (unsigned int j = 0; j < numStepsEntry; ++ j) {
        if (j == 0) {
            zvalsEntry[j] = -entranceParameter3_m - stepSizeEntry;
        } else {
            zvalsEntry[j] = zvalsEntry[j - 1] + stepSizeEntry;
        }
        calcEngeFunction(zvalsEntry[j] / gap_m,
                         engeCoeffsEntry_m,
                         polyOrderEntry_m,
                         fieldValuesEntry[0][j],
                         fieldValuesEntry[1][j],
                         fieldValuesEntry[2][j]);
        fieldValuesEntry[2][j] *= fieldValuesEntry[0][j];
    }

    for (unsigned int j = 0; j < numStepsExit; ++ j) {
        if (j == 0) {
            zvalsExit[j] = exitParameter1_m - stepSizeExit;
        } else {
            zvalsExit[j] = zvalsExit[j - 1] + stepSizeExit;
        }
        calcEngeFunction(zvalsExit[j] / gap_m,
                         engeCoeffsExit_m,
                         polyOrderExit_m,
                         fieldValuesExit[0][j],
                         fieldValuesExit[1][j],
                         fieldValuesExit[2][j]);
        fieldValuesExit[2][j] *= fieldValuesExit[0][j];
    }

    for (unsigned int i = 0; i < 3u; ++ i) {
        gsl_spline_init(entryFieldValues_m[i], &(zvalsEntry[0]), &(fieldValuesEntry[i][0]), numStepsEntry);
        gsl_spline_init(exitFieldValues_m[i], &(zvalsExit[0]), &(fieldValuesExit[i][0]), numStepsExit);
    }
}

void Bend2D::setBendEffectiveLength(double startField, double endField) {

    // Find initial angle.
    double actualBendAngle = calculateBendAngle();

    // Adjust field map to match bend angle.
    double error = std::abs(actualBendAngle - angle_m);
    if(error > 1.0e-6)
        findBendEffectiveLength(startField, endField);

}

void Bend2D::setBendStrength() {

    // Estimate bend field magnitude.
    double mass = RefPartBunch_m->getM();
    double gamma = designEnergy_m / mass + 1.0;
    double betaGamma = sqrt(pow(gamma, 2.0) - 1.0);
    double charge = RefPartBunch_m->getQ();

    fieldAmplitude_m = ((charge / std::abs(charge)) * betaGamma * mass /
                        (Physics::c * designRadius_m));


    // Search for angle if initial guess is not good enough.
    findBendStrength(mass, gamma, betaGamma, charge);
}

void Bend2D::setEngeOriginDelta(double delta) {
    /*
     * This function is used to shift the perpendicular distance of the
     * entrance and exit Enge function origins with respect to the entrance
     * and exit points in the magnet. A positive delta shifts them towards
     * the center of the magnet.
     */
    findChordLength(chordLength_m);

    entranceParameter1_m = delta - std::abs(entranceParameter1_m
                                            - entranceParameter2_m);
    entranceParameter3_m = delta + std::abs(entranceParameter2_m
                                            - entranceParameter3_m);
    entranceParameter2_m = delta;

    exitParameter1_m = -delta - std::abs(exitParameter1_m - exitParameter2_m);
    exitParameter3_m = -delta + std::abs(exitParameter2_m - exitParameter3_m);
    exitParameter2_m = -delta;

    setupFringeWidths();
}

void Bend2D::setFieldCalcParam() {

    cosEntranceAngle_m = cos(entranceAngle_m);
    sinEntranceAngle_m = sin(entranceAngle_m);
    tanEntranceAngle_m = tan(entranceAngle_m);

    deltaBeginEntry_m = std::abs(entranceParameter1_m - entranceParameter2_m);
    deltaEndEntry_m   = std::abs(entranceParameter2_m - entranceParameter3_m);

    deltaBeginExit_m = std::abs(exitParameter1_m - exitParameter2_m);
    deltaEndExit_m   = std::abs(exitParameter2_m - exitParameter3_m);

    double bendAngle = getBendAngle();
    double entranceAngle = getEntranceAngle();

    double rotationAngleAboutZ = getRotationAboutZ();
    Quaternion_t rotationAboutZ(cos(0.5 * rotationAngleAboutZ),
                                sin(0.5 * rotationAngleAboutZ) * Vector_t(0, 0, 1));

    Vector_t rotationAxis(0, -1, 0);
    Quaternion_t halfRotationAboutAxis(cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                                       sin(0.5 * (0.5 * bendAngle - entranceAngle)) * rotationAxis);
    Quaternion_t exitFaceRotation(cos(0.5 * (bendAngle - entranceAngle - exitAngle_m)),
                                  sin(0.5 * (bendAngle - entranceAngle - exitAngle_m)) * rotationAxis);
    Vector_t chord = getChordLength() * halfRotationAboutAxis.rotate(Vector_t(0, 0, 1));
    beginToEnd_lcs_m = CoordinateSystemTrafo(chord, exitFaceRotation.conjugate());
    beginToEnd_m = beginToEnd_lcs_m * CoordinateSystemTrafo(Vector_t(0.0), rotationAboutZ.conjugate());
    toEntranceRegion_m = CoordinateSystemTrafo(Vector_t(0, 0, entranceParameter2_m),
                                               Quaternion(0, 0, 1, 0));
    const CoordinateSystemTrafo fromEndToExitRegion(Vector_t(0, 0, exitParameter2_m),
                                                    Quaternion(1, 0, 0, 0));
    toExitRegion_m = CoordinateSystemTrafo(fromEndToExitRegion *
                                           getBeginToEnd_local());

    Vector_t rotationCenter = Vector_t(-designRadius_m * cosEntranceAngle_m, 0.0, designRadius_m * sinEntranceAngle_m);

    Vector_t maxAngleEntranceAperture;
    if (rotationCenter(2) < 0.0) {
        double tau, tmp;
        Vector_t P(0.5 * aperture_m.second[0], 0.0, 0.0), &R = rotationCenter;
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m + 0.5 * aperture_m.second[0], 2),
                                 &tau,
                                 &tmp);
        tau = (std::abs(1.0 - tmp) < std::abs(1.0 - tau)? tmp: tau);
        maxAngleEntranceAperture = tau * P;
    } else {
        double tau, tmp;
        Vector_t P(-0.5 * aperture_m.second[0], 0.0, 0.0), &R = rotationCenter;
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m - 0.5 * aperture_m.second[0], 2),
                                 &tau,
                                 &tmp);
        tau = (std::abs(1.0 - tmp) < std::abs(1.0 - tau)? tmp: tau);
        maxAngleEntranceAperture = tau * P;
    }
    Vector_t maxAngleExitAperture;
    if (getBeginToEnd_local().transformTo(rotationCenter)(2) > 0.0) {
        double tau, tmp;
        Vector_t P(0.5 * aperture_m.second[0], 0.0, 0.0), R = getBeginToEnd_local().transformTo(rotationCenter);
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m + 0.5 * aperture_m.second[0], 2),
                                 &tau,
                                 &tmp);
        tau = (std::abs(1.0 - tmp) < std::abs(1.0 - tau)? tmp: tau);
        maxAngleExitAperture = getBeginToEnd_local().transformFrom(tau * P);
    } else {
        double tau, tmp;
        Vector_t P(-0.5 * aperture_m.second[0], 0.0, 0.0), R = getBeginToEnd_local().transformTo(rotationCenter);
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m - 0.5 * aperture_m.second[0], 2),
                                 &tau,
                                 &tmp);
        tau = (std::abs(1.0 - tmp) < std::abs(1.0 - tau)? tmp: tau);
        maxAngleExitAperture = getBeginToEnd_local().transformFrom(tau * P);
    }

    maxAngleEntranceAperture -= rotationCenter;

    Quaternion rotation = getQuaternion(maxAngleEntranceAperture, Vector_t(0, 0, 1));
    computeAngleTrafo_m = CoordinateSystemTrafo(rotationCenter,
                                                rotation);
    Vector_t tmp = computeAngleTrafo_m.transformTo(maxAngleExitAperture);
    maxAngle_m = fmod(Physics::two_pi - atan2(tmp(0), tmp(2)), Physics::two_pi);

    maxAngleExitAperture -= rotationCenter;
}

void Bend2D::setGapFromFieldMap() {

    if(gap_m <= 0.0)
        gap_m = fieldmap_m->getFieldGap();
    else if(gap_m != fieldmap_m->getFieldGap())
        adjustFringeFields(gap_m / fieldmap_m->getFieldGap());

}

bool Bend2D::setupBendGeometry(Inform &msg, double &startField, double &endField) {

    chordLength_m = 0.0;
    if(!findChordLength(chordLength_m))
        return false;

    if(treatAsDrift(msg, chordLength_m)) {
        startField_m = startField;
        endField_m = startField + chordLength_m;
        return true;
    }

    reinitialize_m = findIdealBendParameters(chordLength_m);
    if (getType() == RBEND) {
        setEntranceAngle(getEntranceAngle());
    }

    /*
     * Set field map geometry.
     */
    if(aperture_m.second[0] <= 0.0)
        aperture_m.second[0] = designRadius_m / 2.0;
    setFieldCalcParam();

    /*
     * If we are using the default field map, then the origins for the
     * Enge functions will be shifted so that we get the desired bend
     * angle for the given bend strength. (We match the effective length
     * of our field map to the ideal bend length.)
     *
     * If we are not using the default field map, we assume it cannot
     * change so we either leave everything alone (if the user defines
     * the bend strength) or we adjust the bend field to get the right
     * angle.
     */
    elementEdge_m = startField;
    setFieldBoundaries(startField, endField);

    if(fileName_m != "1DPROFILE1-DEFAULT") {
        if(reinitialize_m)
            setBendStrength();
    } else {
        setBendEffectiveLength(startField, endField);
    }

    startField = startField_m;
    endField = endField_m;
    return true;

}

bool Bend2D::setupDefaultFieldMap(Inform &/*msg*/) {

    if(length_m <= 0.0) {
        ERRORMSG("If using \"1DPROFILE1-DEFAULT\" field map you must set the "
                 "chord length of the bend using the L attribute in the OPAL "
                 "input file."
                 << endl);
        return false;
    } else {
        // msg << level2;
        // fieldmap_m->getInfo(&msg);
        return true;
    }

}

void Bend2D::setFieldBoundaries(double startField, double /*endField*/) {

    startField_m = startField - deltaBeginEntry_m / cos(entranceAngle_m);
    endField_m = (startField + angle_m * designRadius_m +
                  deltaEndExit_m / cos(exitAngle_m));

}

void Bend2D::setupPusher(PartBunchBase<double, 3> *bunch) {

    RefPartBunch_m = bunch;
    pusher_m.initialise(bunch->getReference());

}

bool Bend2D::treatAsDrift(Inform &/*msg*/, double chordLength) {
    if(designEnergy_m <= 0.0) {
        WARNMSG("Warning: bend design energy is zero. Treating as drift."
                << endl);
        designRadius_m = 0.0;
        return true;
    } else if(angle_m == 0.0) {

        double refMass = RefPartBunch_m->getM();
        double refGamma = designEnergy_m / refMass + 1.0;
        double refBetaGamma = sqrt(pow(refGamma, 2.0) - 1.0);

        double amplitude = std::abs(fieldAmplitude_m);
        double radius = std::abs(refBetaGamma * refMass / (Physics::c * amplitude));
        double sinArgument = chordLength / (2.0 * radius);

        if(std::abs(sinArgument) > 1.0) {
            WARNMSG("Warning: bend strength and defined reference trajectory "
                    << "chord length are not consistent. Treating bend as drift."
                    << endl);
            designRadius_m = 0.0;
            return true;
        } else
            return false;

    } else if(angle_m == 0.0 &&
              pow(bX_m, 2.0) + pow(bY_m, 2.0) == 0.0) {

        WARNMSG("Warning bend angle/strength is zero. Treating as drift."
                << endl);
        designRadius_m = 0.0;
        return true;

    } else
        return false;
}

std::vector<Vector_t> Bend2D::getOutline() const {
    std::vector<Vector_t> outline;
    Vector_t rotationCenter = Vector_t(-designRadius_m * cosEntranceAngle_m, 0.0, designRadius_m * sinEntranceAngle_m);
    unsigned int numSteps = 2;

    outline.push_back(Vector_t(-0.5 * widthEntranceFringe_m, 0.0, 0.0));
    outline.push_back(Vector_t(-0.5 * widthEntranceFringe_m + entranceParameter1_m * tanEntranceAngle_m, 0.0, entranceParameter1_m));
    outline.push_back(Vector_t(entranceParameter1_m * tanEntranceAngle_m, 0.0, entranceParameter1_m));
    outline.push_back(Vector_t(0.5 * widthEntranceFringe_m + entranceParameter1_m * tanEntranceAngle_m, 0.0, entranceParameter1_m));
    outline.push_back(Vector_t(0.5 * widthEntranceFringe_m, 0.0, 0.0));

    {
        double tau1, tau2;
        Vector_t P(0.5 * aperture_m.second[0], 0.0, 0.0), R = rotationCenter;
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m + 0.5 * aperture_m.second[0], 2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t upperCornerAtEntry = tau1 * P;

        R = getBeginToEnd_local().transformTo(rotationCenter);
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m + 0.5 * aperture_m.second[0], 2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t upperCornerAtExit = getBeginToEnd_local().transformFrom(tau1 * P);

        Quaternion rotation = getQuaternion(upperCornerAtEntry - rotationCenter, Vector_t(0,0,1));
        Vector_t tmp = CoordinateSystemTrafo(rotationCenter, rotation).transformTo(upperCornerAtExit);
        double totalAngle = -fmod(Physics::two_pi - atan2(tmp(0), tmp(2)), Physics::two_pi);
        numSteps = std::max(2.0, std::ceil(-totalAngle / (5.0 * Physics::deg2rad)));
        double dAngle = 0.5 * totalAngle / (1.0 * numSteps - 1.0);

        outline.push_back(upperCornerAtEntry);
        double angle = 0.0;
        for (unsigned int i = 0; i < numSteps; ++ i, angle += dAngle) {
            Quaternion rot(cos(angle), 0.0, sin(angle), 0.0);
            outline.push_back(rot.rotate(upperCornerAtEntry - rotationCenter) + rotationCenter);
        }
        outline.push_back(upperCornerAtExit);
    }

    outline.push_back(getBeginToEnd_local().transformFrom(Vector_t(0.5 * widthExitFringe_m, 0.0, 0.0)));
    outline.push_back(getBeginToEnd_local().transformFrom(Vector_t(0.5 * widthExitFringe_m - exitParameter3_m * tanExitAngle_m, 0.0, exitParameter3_m)));
    outline.push_back(getBeginToEnd_local().transformFrom(Vector_t(-exitParameter3_m * tanExitAngle_m, 0.0, exitParameter3_m)));
    outline.push_back(getBeginToEnd_local().transformFrom(Vector_t(-0.5 * widthExitFringe_m - exitParameter3_m * tanExitAngle_m, 0.0, exitParameter3_m)));
    outline.push_back(getBeginToEnd_local().transformFrom(Vector_t(-0.5 * widthExitFringe_m, 0.0, 0.0)));

    {
        double tau1, tau2;
        Vector_t P(-0.5 * aperture_m.second[0], 0.0, 0.0), R = rotationCenter;
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m - 0.5 * aperture_m.second[0], 2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t lowerCornerAtEntry = tau1 * P;

        R = getBeginToEnd_local().transformTo(rotationCenter);
        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m - 0.5 * aperture_m.second[0], 2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t lowerCornerAtExit = getBeginToEnd_local().transformFrom(tau1 * P);

        Quaternion rotation = getQuaternion(lowerCornerAtEntry - rotationCenter, Vector_t(0,0,1));
        Vector_t tmp = CoordinateSystemTrafo(rotationCenter, rotation).transformTo(lowerCornerAtExit);
        double totalAngle = -fmod(Physics::two_pi - atan2(tmp(0), tmp(2)), Physics::two_pi);
        double dAngle = 0.5 * totalAngle / (1.0 * numSteps - 1.0);

        outline.push_back(lowerCornerAtExit);
        double angle = 0.5 * totalAngle;
        for (unsigned int i = 0; i < numSteps; ++ i, angle -= dAngle) {
            Quaternion rot(cos(angle), 0.0, sin(angle), 0.0);
            outline.push_back(rot.rotate(lowerCornerAtEntry - rotationCenter) + rotationCenter);
        }
        outline.push_back(lowerCornerAtEntry);
    }

    std::ofstream outlineOutput;
    if (Options::writeBendTrajectories && Ippl::myNode() == 0) {
        outlineOutput.open("data/" + OpalData::getInstance()->getInputBasename() + "_" + getName() + "_outline.dat");
        outlineOutput.precision(8);

        for (auto a: outline) {
            outlineOutput << std::setw(18) << a(2)
                          << std::setw(18) << a(0)
                          << "\n";
        }
        outlineOutput << std::setw(18) << outline.front()(2)
                      << std::setw(18) << outline.front()(0)
                      << std::endl;
    }
    return outline;
}

MeshData Bend2D::getSurfaceMesh() const {
    MeshData mesh;
    const Vector_t hgap(0, 0.5 * getFullGap(), 0);
    std::vector<Vector_t> outline = getOutline();

    unsigned int size = outline.size();
    unsigned int last = size - 1;
    unsigned int numSteps = (size - 14) / 2;
    unsigned int midIdx = numSteps + 7;

    for (unsigned int i = 0; i < 6; ++ i) {
        mesh.vertices_m.push_back(outline[i] + 2 * hgap);
    }
    for (unsigned int i = 6; i < 6 + numSteps; ++ i) {
        mesh.vertices_m.push_back(outline[i] + hgap);
    }
    for (unsigned int i = 6 + numSteps; i < 13 + numSteps; ++ i) {
        mesh.vertices_m.push_back(outline[i] + 2 * hgap);
    }
    for (unsigned int i = 13 + numSteps; i < 13 + 2 * numSteps; ++ i) {
        mesh.vertices_m.push_back(outline[i] + hgap);
    }
    mesh.vertices_m.push_back(outline[last] + 2 * hgap);

    for (unsigned int i = 0; i < 6; ++ i) {
        mesh.vertices_m.push_back(outline[i] - 2 * hgap);
    }
    for (unsigned int i = 6; i < 6 + numSteps; ++ i) {
        mesh.vertices_m.push_back(outline[i] - hgap);
    }
    for (unsigned int i = 6 + numSteps; i < 13 + numSteps; ++ i) {
        mesh.vertices_m.push_back(outline[i] - 2 * hgap);
    }
    for (unsigned int i = 13 + numSteps; i < 13 + 2 * numSteps; ++ i) {
        mesh.vertices_m.push_back(outline[i] - hgap);
    }
    mesh.vertices_m.push_back(outline[last] - 2 * hgap);

    mesh.vertices_m.push_back(outline[0]);
    mesh.vertices_m.push_back(outline[4]);
    mesh.vertices_m.push_back(outline[midIdx]);
    mesh.vertices_m.push_back(outline[midIdx + 4]);

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 1, 0));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 4, 3));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 5, 4));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 0, last));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, last, 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(last, last - 1, 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(5, last - 1, 6));

    // exit region top
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx + 2, midIdx + 1, midIdx));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx + 2, midIdx + 4, midIdx + 3));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx + 2, midIdx + 5, midIdx + 4));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx + 2, midIdx, midIdx - 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx + 2, midIdx - 1, midIdx + 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx - 1, midIdx - 2, midIdx + 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(midIdx + 5, midIdx - 2, midIdx + 6));

    // entry region bottom
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 2, size + 0, size + 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 2, size + 3, size + 4));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 2, size + 4, size + 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 2, size + last, size + 0));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 2, size + 5, size + last));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + last, size + 5, size + last - 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 5, size + 6, size + last - 1));

    // exit region bottom
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx + 2, size + midIdx + 0, size + midIdx + 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx + 2, size + midIdx + 3, size + midIdx + 4));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx + 2, size + midIdx + 4, size + midIdx + 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx + 2, size + midIdx - 1, size + midIdx + 0));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx + 2, size + midIdx + 5, size + midIdx - 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx - 1, size + midIdx + 5, size + midIdx - 2));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + midIdx + 5, size + midIdx + 6, size + midIdx - 2));

    // central region
    for (unsigned int i = 6; i < 5 + numSteps; ++ i) {
        mesh.triangles_m.push_back(Vektor<unsigned int, 3>(i, last + 5 - i, i + 1));
        mesh.triangles_m.push_back(Vektor<unsigned int, 3>(last + 5 - i, last + 4 - i, i + 1));

        mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + i, size + i + 1, size + last + 5 - i));
        mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + last + 5 - i, size + i + 1, size + last + 4 - i));
    }

    // side
    for (unsigned int i = 0; i < size - 2; ++ i) {
        if (i == 4u || i == 5u ||
            i == 5 + numSteps || i == 6 + numSteps ||
            i == 11 + numSteps || i == 12 + numSteps) continue;

        unsigned int next = (i + 1) % size;
        mesh.triangles_m.push_back(Vektor<unsigned int, 3>(i, next, next + size));
        mesh.triangles_m.push_back(Vektor<unsigned int, 3>(i, next + size, i + size));
    }

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(last, 0, last-1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2*size, last - 1, 0));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2*size, size + last - 1, last - 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2*size, size, size + last - 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + last, size + last - 1, size));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 5, 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 6, 2*size + 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2*size + 1, 6, size + 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 4, 2*size + 1, size + 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 4, size + 6, size + 5));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(6 + numSteps, midIdx, 5 + numSteps));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(5 + numSteps, midIdx, 2*size + 2));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(5 + numSteps, 2*size + 2, size + 5 + numSteps));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 5 + numSteps, 2*size + 2, size + midIdx));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 6 + numSteps, size + 5 + numSteps, size + midIdx));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(6 + midIdx, 4 + midIdx, 5 + midIdx));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(6 + midIdx, 2*size + 3, 4 + midIdx));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2*size + 3, 6 + midIdx, size + 6 + midIdx));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 4 + midIdx,  2*size + 3, size + 6 + midIdx));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(size + 4 + midIdx, size + 6 + midIdx, size + 5 + midIdx));

    Vector_t rotationCenter(-designRadius_m * cosEntranceAngle_m, 0.0, designRadius_m * sinEntranceAngle_m);

    Vector_t P1 = toEntranceRegion_m.transformFrom(Vector_t(0, 0, entranceParameter1_m));
    Vector_t P2 = toExitRegion_m.transformFrom(Vector_t(0, 0, exitParameter1_m));

    Vector_t T = cross(P1 - rotationCenter, P2 - rotationCenter);
    if (T[1] > 0) { // fringe fields are overlapping
        Vector_t dir1 = toEntranceRegion_m.rotateFrom(Vector_t(1, 0, 0));
        if (this->getType() == ElementBase::RBEND ||
            std::abs(entranceAngle_m + exitAngle_m - angle_m) < 1e-8) {
            mesh.decorations_m.push_back(std::make_pair(0.5 * (P1 + P2) - 0.25 * dir1,
                                                        0.5 * (P1 + P2) + 0.25 * dir1));
        } else {
            Vector_t dir2 = toExitRegion_m.rotateFrom(Vector_t(-1, 0, 0));
            Tenzor<double, 3> inv;
            double det = -dir1[0] * dir2[2] + dir1[2] * dir2[0];
            inv(0, 0) = -dir2[2] / det;
            inv(0, 2) = dir2[0] / det;
            inv(1,1) = 1.0;
            inv(2, 0) = -dir1[2] / det;
            inv(2, 2) = dir1[0] / det;
            Vector_t Tau = dot(inv, P2 - P1);
            Vector_t crossPoint = P1 + Tau[0] * dir1;
            double angle = asin(cross(dir1, dir2)[1]);
            Quaternion halfRot(cos(0.25 * angle), sin(0.25 * angle) * Vector_t(0, 1, 0));
            Vector_t P = halfRot.rotate(dir1);
            Vector_t R = crossPoint - rotationCenter;

            double radius = designRadius_m - 0.5 * aperture_m.second[0];
            double tau1, tau2;
            gsl_poly_solve_quadratic(dot(P,P),
                                     2 * dot(P,R),
                                     dot(R,R) - std::pow(radius, 2),
                                     &tau1,
                                     &tau2);
            if (euclidean_norm(crossPoint + tau1 * P - P1) > euclidean_norm(crossPoint + tau2 * P - P1)) {
                std::swap(tau1, tau2);
            }
            Vector_t lowerCornerFringeLimit = crossPoint + tau1 * P;

            radius = designRadius_m + 0.5 * aperture_m.second[0];
            gsl_poly_solve_quadratic(dot(P,P),
                                     2 * dot(P,R),
                                     dot(R,R) - std::pow(radius, 2),
                                     &tau1,
                                     &tau2);
            if (euclidean_norm(crossPoint + tau1 * P - P1) > euclidean_norm(crossPoint + tau2 * P - P1)) {
                std::swap(tau1, tau2);
            }
            Vector_t upperCornerFringeLimit = crossPoint + tau1 * P;

            mesh.decorations_m.push_back(std::make_pair(lowerCornerFringeLimit,
                                                        upperCornerFringeLimit));
        }
    } else {

        double tau1, tau2;
        Vector_t P(-0.5 * aperture_m.second[0], 0.0, 0.0);
        Vector_t R = Vector_t(0, 0, entranceParameter3_m) - rotationCenter;
        gsl_poly_solve_quadratic(dot(P,P),
                                 2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m - 0.5 * aperture_m.second[0],2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t lowerCornerFringeLimitEntrance = Vector_t(0, 0, entranceParameter3_m) + tau1 * P;

        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m + 0.5 * aperture_m.second[0],2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t upperCornerFringeLimitEntrance = Vector_t(0, 0, entranceParameter3_m) - tau1 * P;

        P = Vector_t(0.5 * aperture_m.second[0], 0.0, 0.0);
        R = Vector_t(0, 0, exitParameter1_m) - getBeginToEnd_local().transformTo(rotationCenter);
        gsl_poly_solve_quadratic(dot(P,P),
                                 2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m + 0.5 * aperture_m.second[0],2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t upperCornerFringeLimitExit = getBeginToEnd_local().transformFrom(Vector_t(0, 0, exitParameter1_m) + tau1 * P);

        gsl_poly_solve_quadratic(dot(P,P),
                                 -2 * dot(P,R),
                                 dot(R,R) - std::pow(designRadius_m - 0.5 * aperture_m.second[0],2),
                                 &tau1,
                                 &tau2);
        tau1 = (std::abs(1.0 - tau2) < std::abs(1.0 - tau1)? tau2: tau1);
        Vector_t lowerCornerFringeLimitExit = getBeginToEnd_local().transformFrom(Vector_t(0, 0, exitParameter1_m) - tau1 * P);

        mesh.decorations_m.push_back(std::make_pair(lowerCornerFringeLimitEntrance, upperCornerFringeLimitEntrance));
        mesh.decorations_m.push_back(std::make_pair(lowerCornerFringeLimitExit, upperCornerFringeLimitExit));

    }

    mesh.decorations_m.push_back(std::make_pair(Vector_t(entranceParameter1_m * tanEntranceAngle_m, 0.0, entranceParameter1_m),
                                                Vector_t(0.0)));
    mesh.decorations_m.push_back(std::make_pair(getBeginToEnd_local().transformFrom(Vector_t(0.0)),
                                                getBeginToEnd_local().transformFrom(Vector_t(-exitParameter3_m * tanExitAngle_m, 0.0, exitParameter3_m))));

    return mesh;
}

bool Bend2D::isInside(const Vector_t &R) const {
    if (std::abs(R(1)) > gap_m) return false;

    if (inMagnetEntranceRegion(R)) {
        return true;
    }

    if (inMagnetExitRegion(R)) {
        return true;
    }

    return (std::abs(R(1)) < 0.5 * gap_m && inMagnetCentralRegion(R));
}

void Bend2D::setupFringeWidths()
{
    widthEntranceFringe_m = 2 * std::min(entranceParameter3_m - entranceParameter1_m, aperture_m.second[0]) + aperture_m.second[0];
    widthExitFringe_m = 2 * std::min(exitParameter3_m - exitParameter1_m, aperture_m.second[0]) + aperture_m.second[0];
}

//set the number of slices for map tracking
void Bend2D::setNSlices(const std::size_t& nSlices) {
    nSlices_m = nSlices;
}

//get the number of slices for map tracking
std::size_t Bend2D::getNSlices() const {
    return nSlices_m;
}

/// Get entrance fringe field length.
//  Used to create fringe fields in ThickTracker, (before edge[m], after edge[m])
std::array<double,2> Bend2D::getEntranceFringeFieldLength() const {
    double entrParam1, entrParam2, entrParam3;

    fieldmap_m->get1DProfile1EntranceParam(entrParam1,
                entrParam2,
                entrParam3);

    std::array<double,2> entFFL;
    entFFL[0]=entrParam2-entrParam1; //before edge
    entFFL[1]=entrParam3-entrParam2; //after edge


    entFFL[0] = ( entranceParameter2_m - entranceParameter1_m ); //before edge
    entFFL[1] = ( entranceParameter3_m - entranceParameter2_m ); //after edge

    return entFFL;
}

/// Get exit fringe field length.
//  Used to create fringe fields in ThickTracker, (after edge[m], before edge[m])
std::array<double,2> Bend2D::getExitFringeFieldLength() const {
    std::array<double,2> extFFL;

    double exitParam1, exitParam2, exitParam3;

    fieldmap_m->get1DProfile1ExitParam(exitParam1,
                exitParam2,
                exitParam3);

    extFFL[0]=exitParam3-exitParam2; //after edge
    extFFL[1]=exitParam2-exitParam1; //before edge

    extFFL[0] = ( exitParameter3_m-exitParameter2_m ); //after edge
    extFFL[1] = ( exitParameter2_m-exitParameter1_m ); //before edge

    return extFFL;
}

