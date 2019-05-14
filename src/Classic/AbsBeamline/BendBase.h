#ifndef CLASSIC_BENDBASE_H
#define CLASSIC_BENDBASE_H

#include "AbsBeamline/Component.h"

#include <vector>
#include <string>

class BendBase: public Component {
public:
    BendBase();
    BendBase(const BendBase &);
    BendBase(const std::string &);

    /// Indicates that element bends the beam.
    virtual bool bends() const;

    void setLength(double length);
    double getLength() const;
    double getChordLength() const;
    virtual void setBendAngle(double angle);
    double getBendAngle() const;
    virtual void setEntranceAngle(double entranceAngle);
    double getEntranceAngle() const;
    virtual double getExitAngle() const = 0;
    void setFullGap(double);
    double getFullGap() const;

    virtual void setDesignEnergy(const double& energy, bool changeable = true);
    double getDesignEnergy() const;
    std::vector<Vector_t> getDesignPath() const;

    void setFieldAmplitude(double k0, double k0s);
    double getFieldAmplitude() const;

    void setFieldMapFN(std::string fileName);
    std::string getFieldMapFN() const;

protected:
    double length_m;
    double chordLength_m;
    double angle_m;
    double entranceAngle_m;     /// Angle between incoming reference trajectory

    double gap_m;

    double designEnergy_m;      /// Bend design energy (eV).
    bool designEnergyChangeable_m;
    /// Map of reference particle trajectory.
    std::vector<Vector_t> refTrajMap_m;

    double bX_m;
    double bY_m;
    double fieldAmplitude_m;

    std::string fileName_m;
};

inline
bool BendBase::bends() const {
    return true;
}

inline
void BendBase::setLength(double length) {
    length_m = std::abs(length);
}

inline
double BendBase::getLength() const
{
    return length_m;
}

inline
double BendBase::getChordLength() const {
    return chordLength_m;
}

inline
void BendBase::setBendAngle(double angle) {
    angle_m = angle;
}

inline
double BendBase::getBendAngle() const {
    return angle_m;
}

inline
void BendBase::setEntranceAngle(double angle)
{
    entranceAngle_m = angle;
}

inline
double BendBase::getEntranceAngle() const {
    return entranceAngle_m;
}

inline
void BendBase::setFullGap(double gap) {
    gap_m = std::abs(gap);
}

inline
double BendBase::getFullGap() const {
    return gap_m;
}

inline
void BendBase::setDesignEnergy(const double& energy, bool changeable) {
    if (designEnergyChangeable_m) {
        designEnergy_m = std::abs(energy) * 1e6;
        designEnergyChangeable_m = changeable;
    }
}

inline
double BendBase::getDesignEnergy() const {
    return designEnergy_m;
}

inline
void BendBase::setFieldAmplitude(double k0, double k0s) {
    bY_m = k0;
    bX_m = k0s;
}

inline
double BendBase::getFieldAmplitude() const
{
    return fieldAmplitude_m;
}

inline
void BendBase::setFieldMapFN(std::string fileName) {
    fileName_m = fileName;
}

inline
std::string BendBase::getFieldMapFN() const {
    return fileName_m;
}


#endif