#include "AbsBeamline/BendBase.h"

#include "Algorithms/PartBunchBase.h"

#include <cmath>

BendBase::BendBase():
    BendBase("")
{}

BendBase::BendBase(const BendBase &right):
    Component(right),
    length_m(right.length_m),
    chordLength_m(right.chordLength_m),
    angle_m(right.angle_m),
    entranceAngle_m(right.entranceAngle_m),
    fieldmap_m(right.fieldmap_m),
    gap_m(right.gap_m),
    designEnergy_m(right.designEnergy_m),
    designEnergyChangeable_m(true),
    refTrajMap_m(right.refTrajMap_m),
    fieldAmplitudeX_m(right.fieldAmplitudeX_m),
    fieldAmplitudeY_m(right.fieldAmplitudeY_m),
    fieldAmplitude_m(right.fieldAmplitude_m),
    fileName_m(right.fileName_m)
{}

BendBase::BendBase(const std::string &name):
    Component(name),
    length_m(0.0),
    chordLength_m(0.0),
    angle_m(0.0),
    entranceAngle_m(0.0),
    fieldmap_m(NULL),
    gap_m(0.0),
    designEnergy_m(0.0),
    designEnergyChangeable_m(true),
    fieldAmplitudeX_m(0.0),
    fieldAmplitudeY_m(0.0),
    fieldAmplitude_m(0.0),
    fileName_m("")
{}


std::vector<Vector_t> BendBase::getDesignPath() const {
    unsigned int size = refTrajMap_m.size();
    std::vector<Vector_t> designPath(size);
    // double angleZ = getRotationAboutZ();
    // Quaternion rotationAboutZ(cos(angleZ / 2), sin(angleZ / 2) * Vector_t(0, 0, 1));
    for (unsigned int i = 0; i < size; ++ i) {
        Vector_t currentPosition = refTrajMap_m[i];
        designPath[i] = currentPosition;//rotationAboutZ.rotate(currentPosition);
    }

    return designPath;
}

void BendBase::setFieldAmplitude(double k0, double k0s) {
    fieldAmplitudeY_m = k0;
    fieldAmplitudeX_m = k0s;
}

double BendBase::calcDesignRadius(double fieldAmplitude) const
{
    double mass = RefPartBunch_m->getM();
    double betaGamma = calcBetaGamma();
    double charge = RefPartBunch_m->getQ();
    // Lorentz force: condition for a circular orbit
    return std::abs(betaGamma * mass / (Physics::c * fieldAmplitude * charge));
}

double BendBase::calcFieldAmplitude(double radius) const
{
    double mass = RefPartBunch_m->getM();
    double betaGamma = calcBetaGamma();
    double charge = RefPartBunch_m->getQ();
    // Lorentz force: condition for a circular orbit
    return betaGamma * mass / (Physics::c * radius * charge);
}

double BendBase::calcBendAngle(double chordLength, double radius) const
{
    return 2.0 * std::asin(chordLength / (2.0 * radius));
}

double BendBase::calcDesignRadius(double chordLength, double angle) const
{
    return chordLength / (2.0 * std::sin(angle / 2.0));
}

double BendBase::calcGamma() const
{
    double mass = RefPartBunch_m->getM();
    return designEnergy_m / mass + 1.0;
}

double BendBase::calcBetaGamma() const
{
    double gamma = calcGamma();
    return std::sqrt(std::pow(gamma, 2.0) - 1.0);
}