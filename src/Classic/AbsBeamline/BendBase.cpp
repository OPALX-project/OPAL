#include "AbsBeamline/BendBase.h"

BendBase::BendBase():
    BendBase("")
{}

BendBase::BendBase(const BendBase &right):
    Component(right),
    length_m(right.length_m),
    chordLength_m(right.chordLength_m),
    angle_m(right.angle_m),
    entranceAngle_m(right.entranceAngle_m),
    gap_m(right.gap_m),
    designEnergy_m(right.designEnergy_m),
    designEnergyChangeable_m(true),
    refTrajMap_m(right.refTrajMap_m),
    bX_m(right.bX_m),
    bY_m(right.bY_m),
    fieldAmplitude_m(right.fieldAmplitude_m),
    fileName_m(right.fileName_m)
{
    setElType(isDipole);
}

BendBase::BendBase(const std::string &name):
    Component(name),
    length_m(0.0),
    chordLength_m(0.0),
    angle_m(0.0),
    entranceAngle_m(0.0),
    gap_m(0.0),
    designEnergy_m(0.0),
    designEnergyChangeable_m(true),
    bX_m(0.0),
    bY_m(0.0),
    fieldAmplitude_m(0.0),
    fileName_m("")
{
    setElType(isDipole);
}


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