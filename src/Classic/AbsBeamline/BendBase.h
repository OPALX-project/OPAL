#ifndef CLASSIC_BENDBASE_H
#define CLASSIC_BENDBASE_H

#include "AbsBeamline/Component.h"

#include <vector>
#include <string>

class Fieldmap;

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
    /// Calculate design radius from design energy and field amplitude
    double calcDesignRadius(double fieldAmplitude) const;
    /// Calculate field amplitude from design energy and radius
    double calcFieldAmplitude(double radius) const;
    /// Calculate bend angle from chord length and design radius
    double calcBendAngle(double chordLength, double radius) const;
    /// Calculate design radius from chord length and bend angle
    double calcDesignRadius(double chordLength, double angle) const;
    /// Calculate gamma from design energy
    double calcGamma() const;
    /// Calculate beta*gamma from design energy
    double calcBetaGamma() const;

    double length_m;
    double chordLength_m;
    double angle_m;         ///< Bend angle
    double entranceAngle_m; ///< Angle between incoming reference trajectory
                            ///< and the entrance face of the magnet (radians).
    Fieldmap *fieldmap_m;      ///< Magnet field map.
    const bool fast_m = false; ///< Flag to turn on fast field calculation.

    double gap_m; ///< Full vertical gap of the magnets.

    double designEnergy_m; ///< Bend design energy (eV).
    bool designEnergyChangeable_m;
    /// Map of reference particle trajectory.
    std::vector<Vector_t> refTrajMap_m;

    double fieldAmplitudeX_m; ///< Field amplitude in x direction.
                              ///< Value not updated if user defines strength with angle
    double fieldAmplitudeY_m; ///< Field amplitude in y direction.
                              ///< Value not updated if user defines strength with angle

    double fieldAmplitude_m;  ///< Field amplitude.

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