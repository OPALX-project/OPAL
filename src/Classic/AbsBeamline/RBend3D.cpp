// ------------------------------------------------------------------------
// $RCSfile: RBend3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend3D
//   Defines the abstract interface for a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend3D.h"
#include "Algorithms/PartBunchBase.h"
#include "Steppers/BorisPusher.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Structure/MeshGenerator.h"
#include "Physics/Physics.h"
#include "Utilities/Util.h"

#include <iostream>
#include <fstream>
#include <cmath>

extern Inform *gmsg;

// Class RBend3D
// ------------------------------------------------------------------------

RBend3D::RBend3D():
    RBend3D("")
{}


RBend3D::RBend3D(const RBend3D &right):
    BendBase(right),
    fieldAmplitudeError_m(right.fieldAmplitudeError_m),
    startField_m(right.startField_m),
    lengthField_m(right.lengthField_m),
    geometry_m(right.geometry_m),
    dummyField_m() {
}


RBend3D::RBend3D(const std::string &name):
    BendBase(name),
    fieldAmplitudeError_m(0.0),
    startField_m(0.0),
    lengthField_m(0.0),
    geometry_m(),
    dummyField_m() {
}


RBend3D::~RBend3D() {
}

void RBend3D::accept(BeamlineVisitor &visitor) const {
    visitor.visitRBend3D(*this);
}

bool RBend3D::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool RBend3D::apply(const Vector_t &R, const Vector_t &/*P*/, const  double &/*t*/, Vector_t &/*E*/, Vector_t &B) {
    const Vector_t tmpR(R(0), R(1), R(2) - startField_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

    B += (fieldAmplitude_m + fieldAmplitudeError_m) * tmpB;

    return false;
}

bool RBend3D::applyToReferenceParticle(const Vector_t &R, const Vector_t &/*P*/, const  double &/*t*/, Vector_t &/*E*/, Vector_t &B) {
    const Vector_t tmpR(R(0), R(1), R(2) - startField_m);
    Vector_t tmpE(0.0), tmpB(0.0);

    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    B += fieldAmplitude_m * tmpB;

    return false;
}

void RBend3D::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    Inform msg("RBend3D ", *gmsg);

    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(fileName_m, fast_m);
    if(fieldmap_m != NULL) {
        msg << level2 << getName() << " using file ";
        fieldmap_m->getInfo(&msg);
        goOnline(0.0);

        double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
        fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

        if (length_m == 0.0) {
            chordLength_m = 0.0;
            double fieldLength = zEnd - zBegin;
            double z = 0.0, dz = fieldLength / 1000;
            Vector_t E(0.0), B(0.0);
            while (z < fieldLength && B(1) < 0.5) {
                fieldmap_m->getFieldstrength(Vector_t(0.0, 0.0, z), E, B);
                z += dz;
            }
            double zEntryEdge = z;
            z = fieldLength;
            B(1) = 0.0;
            while (z > 0.0 && B(1) < 0.5) {
                fieldmap_m->getFieldstrength(Vector_t(0.0, 0.0, z), E, B);
                z -= dz;
            }
            chordLength_m = z - zEntryEdge;
            length_m = chordLength_m;
        } else {
            chordLength_m = length_m;
        }

        startField_m = zBegin;
        lengthField_m = zEnd - zBegin;
        endField = startField + lengthField_m;

        if (angle_m != 0.0) {
            if (angle_m < 0.0) {
                // Negative angle is a positive bend rotated 180 degrees.
                entranceAngle_m = std::copysign(1, angle_m) * entranceAngle_m;
                angle_m = std::abs(angle_m);
                rotationZAxis_m += Physics::pi;
            }

            Vector_t B(0.0), E(0.0);
            double z = 0.0, dz = lengthField_m / 999;
            double integratedBy = 0.0;

            fieldmap_m->getFieldstrength(Vector_t(0.0, 0.0, z), E, B);
            integratedBy += 0.5 * B(1);
            z = dz;
            while (z < lengthField_m) {
                B = 0.0; E = 0.0;
                fieldmap_m->getFieldstrength(Vector_t(0.0, 0.0, z), E, B);
                integratedBy += B(1);
                z += dz;
            }
            integratedBy -= 0.5 * B(1);
            integratedBy *= lengthField_m / 1000;

            // estimate magnitude of field
            fieldAmplitude_m = calcFieldAmplitude(integratedBy / angle_m);

            double angle = trackRefParticleThrough(bunch->getdT());
            double error = (angle_m - angle) / angle_m;

            for (unsigned int i = 0; i < 10; ++ i) {
                fieldAmplitude_m *= (0.5 * error + 1.0);
                angle = trackRefParticleThrough(bunch->getdT());
                error = (angle_m - angle) / angle_m;

                if (std::abs(error) < 0.001) break;
            }

            if (Options::writeBendTrajectories)
                trackRefParticleThrough(bunch->getdT(), true);
        } else {
            double refCharge = bunch->getQ();
            rotationZAxis_m += std::atan2(fieldAmplitudeX_m, fieldAmplitudeY_m);
            if (refCharge < 0.0) {
                rotationZAxis_m -= Physics::pi;
            }
            fieldAmplitude_m = std::copysign(1.0, refCharge) * std::hypot(fieldAmplitudeX_m, fieldAmplitudeY_m);
            angle_m = trackRefParticleThrough(bunch->getdT(), Options::writeBendTrajectories);
        }

    } else {
        endField = startField;
    }
}

void RBend3D::finalise()
{}

void RBend3D::goOnline(const double &) {
    Fieldmap::readMap(fileName_m);
    online_m = true;
}

void RBend3D::goOffline() {
    Fieldmap::freeMap(fileName_m);
    online_m = false;
}

void RBend3D::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + length_m;
}


ElementBase::ElementType RBend3D::getType() const {
    return RBEND3D;
}

bool RBend3D::isInside(const Vector_t &r) const {
    Vector_t tmpR(r(0), r(1), r(2) - startField_m);
    return fieldmap_m->isInside(tmpR);
}

ElementBase* RBend3D::clone() const {
    return new RBend3D(*this);
}

BGeometryBase& RBend3D::getGeometry() {
    return geometry_m;
}

const BGeometryBase& RBend3D::getGeometry() const {
    return geometry_m;
}

EMField &RBend3D::getField() {
    return dummyField_m;
}

const EMField &RBend3D::getField() const {
    return dummyField_m;
}

MeshData RBend3D::getSurfaceMesh() const {
    Vector_t XIni, XFinal;
    fieldmap_m->getFieldDimensions(XIni(0), XFinal(0),
                                   XIni(1), XFinal(1),
                                   XIni(2), XFinal(2));

    MeshData mesh;
    mesh.vertices_m.push_back(Vector_t(XIni(0),   XIni(1),   XIni(2)));
    mesh.vertices_m.push_back(Vector_t(XIni(0),   XIni(1),   XFinal(2)));
    mesh.vertices_m.push_back(Vector_t(XIni(0),   XFinal(1), XIni(2)));
    mesh.vertices_m.push_back(Vector_t(XIni(0),   XFinal(1), XFinal(2)));
    mesh.vertices_m.push_back(Vector_t(XFinal(0), XIni(1),   XIni(2)));
    mesh.vertices_m.push_back(Vector_t(XFinal(0), XIni(1),   XFinal(2)));
    mesh.vertices_m.push_back(Vector_t(XFinal(0), XFinal(1), XIni(2)));
    mesh.vertices_m.push_back(Vector_t(XFinal(0), XFinal(1), XFinal(2)));

    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(0, 1, 2));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 1, 3));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 6, 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(5, 6, 7));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(1, 5, 3));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(3, 5, 7));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 2, 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(4, 0, 2));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(0, 4, 1));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(1, 4, 5));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(2, 3, 6));
    mesh.triangles_m.push_back(Vektor<unsigned int, 3>(3, 7, 6));

    return mesh;
}

double RBend3D::trackRefParticleThrough(double dt, bool print) {
    const double refGamma     = calcGamma();
    const double refBetaGamma = calcBetaGamma();
    const double stepSize     = refBetaGamma / refGamma * Physics::c * dt;
    const Vector_t scaleFactor(Physics::c * dt);
    print = print && (Ippl::myNode() == 0);

    std::ofstream trajectoryOutput;
    if (print) {
        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            OpalData::getInstance()->getInputBasename() + "_" + getName() + "_traj.dat"
        });
        trajectoryOutput.open(fname);
        trajectoryOutput.precision(12);
        trajectoryOutput << "# " << std::setw(18) << "s"
                         << std::setw(20) << "x"
                         << std::setw(20) << "z"
                         << std::setw(20) << "By"
                         << std::endl;
    }

    double deltaS = 0.0;
    BorisPusher pusher(*RefPartBunch_m->getReference());

    Vector_t X(0.0), P(0.0);
    X(0) = startField_m * std::tan(entranceAngle_m);
    X(2) = startField_m;
    P(0) = refBetaGamma * std::sin(entranceAngle_m);
    P(2) = refBetaGamma * std::cos(entranceAngle_m);

    while ((X(2) - startField_m) < lengthField_m && 0.5 * deltaS < lengthField_m) {
        Vector_t E(0.0), B(0.0);

        X /= scaleFactor;
        pusher.push(X, P, dt);
        X *= scaleFactor;

        applyToReferenceParticle(X, P, 0.0, E, B);
        if (print) {
            trajectoryOutput << std::setw(20) << deltaS + 0.5 * stepSize
                             << std::setw(20) << X(0)
                             << std::setw(20) << X(2)
                             << std::setw(20) << B(1)
                             << std::endl;
        }

        pusher.kick(X, P, E, B, dt);

        X /= scaleFactor;
        pusher.push(X, P, dt);
        X *= scaleFactor;

        deltaS += stepSize;
    }

    return -std::atan2(P(0), P(2)) + entranceAngle_m;
}