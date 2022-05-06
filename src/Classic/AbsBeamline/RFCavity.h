//
// Class RFCavity
//   Defines the abstract interface for for RF cavities.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_RFCavity_HH
#define CLASSIC_RFCavity_HH

#include "AbsBeamline/Component.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Physics/Physics.h"

#include <boost/bimap.hpp>

#include <cmath>
#include <string>

class Fieldmap;

enum class CavityType: unsigned short {
    SW,
    SGSW
};

class RFCavity: public Component {

public:

    /// Constructor with given name.
    explicit RFCavity(const std::string& name);

    RFCavity();
    RFCavity(const RFCavity&);
    virtual ~RFCavity();

    /// Apply visitor to RFCavity.
    virtual void accept(BeamlineVisitor& ) const override;

    /// Get RF amplitude.
    virtual double getAmplitude() const = 0;

    /// Get RF frequency.
    virtual double getFrequency() const = 0;
    void setFrequency(double freq);

    /// Get RF phase.
    virtual double getPhase() const = 0;

    void dropFieldmaps();

    /// Set the name of the field map
    virtual void setFieldMapFN(const std::string& fmapfn);
    virtual std::string getFieldMapFN() const;

    virtual void setAmplitudem(double vPeak);
    virtual double getAmplitudem() const;

    virtual void setAmplitudeError(double vPeakError);
    virtual double getAmplitudeError() const;

    virtual void setFrequencym(double freq);
    virtual double getFrequencym() const;

    virtual void setPhasem(double phase);
    virtual double getPhasem() const;
    double getPhasem(double t) const;

    virtual void setPhaseError(double phaseError);
    virtual double getPhaseError() const;

    void setCavityType(const std::string& type);
    CavityType getCavityType() const;
    std::string getCavityTypeString() const;

    virtual void setFast(bool fast);
    virtual bool getFast() const;

    virtual void setAutophaseVeto(bool veto = true);
    virtual bool getAutophaseVeto() const;

    virtual double getAutoPhaseEstimate(const double& E0, const double& t0, const double& q, const double& m);
    virtual double getAutoPhaseEstimateFallback(double E0, double t0, double q, double m);

    virtual std::pair<double, double> trackOnAxisParticle(const double& p0,
                                                          const double& t0,
                                                          const double& dt,
                                                          const double& q,
                                                          const double& mass,
                                                          std::ofstream *out = nullptr);

    virtual bool apply(const size_t& i,
                       const double& t,
                       Vector_t& E,
                       Vector_t& B) override;

    virtual bool apply(const Vector_t& R,
                       const Vector_t& P,
                       const double& t,
                       Vector_t& E,
                       Vector_t& B) override;

    virtual bool applyToReferenceParticle(const Vector_t& R,
                                          const Vector_t& P,
                                          const double& t,
                                          Vector_t& E,
                                          Vector_t& B) override;

    virtual void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) override;

    virtual void initialise(PartBunchBase<double, 3>* bunch,
                            std::shared_ptr<AbstractTimeDependence> freq_atd,
                            std::shared_ptr<AbstractTimeDependence> ampl_atd,
                            std::shared_ptr<AbstractTimeDependence> phase_atd);

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double& kineticEnergy) override;

    virtual void goOffline() override;

    virtual void setDesignEnergy(const double& ekin, bool changeable = true) override;
    virtual double getDesignEnergy() const override;

    void setRmin(double rmin);
    virtual double getRmin() const;

    void setRmax(double rmax);
    virtual double getRmax() const;

    void setAzimuth(double angle);
    virtual double getAzimuth() const;

    void setPerpenDistance(double pdis);
    virtual double getPerpenDistance() const;

    void setGapWidth(double gapwidth);
    virtual double getGapWidth() const;

    void setPhi0(double phi0);
    virtual double getPhi0() const;

    virtual double getCosAzimuth() const;

    virtual double getSinAzimuth() const;

    virtual double getCycFrequency() const;

    void getMomentaKick(const double normalRadius,
                        double momentum[],
                        const double t,
                        const double dtCorrt,
                        const int PID,
                        const double restMass,
                        const int chargenumber);

    double spline(double z, double *za);

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual bool isInside(const Vector_t& r) const override;

    void setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    void setAmplitudeModelName(std::string name);
    std::string getAmplitudeModelName();

    void setPhaseModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    void setPhaseModelName(std::string name);
    std::string getPhaseModelName();

    void setFrequencyModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    void setFrequencyModelName(std::string name);
    std::string getFrequencyModelName();

    virtual double getElementLength() const override;
    virtual void getElementDimensions(double& begin,
                                      double& end) const override;

    virtual CoordinateSystemTrafo getEdgeToBegin() const override;
    virtual CoordinateSystemTrafo getEdgeToEnd() const override;

protected:
    std::shared_ptr<AbstractTimeDependence> phaseTD_m;
    std::string phaseName_m;
    std::shared_ptr<AbstractTimeDependence> amplitudeTD_m;
    std::string amplitudeName_m;
    std::shared_ptr<AbstractTimeDependence> frequencyTD_m;
    std::string frequencyName_m;

    std::string filename_m;      /**< The name of the inputfile*/

    double scale_m;              /**< scale multiplier*/
    double scaleError_m;         /**< additive scale error*/
    double phase_m;              /**< phase shift of time varying field (rad)*/
    double phaseError_m;         /**< phase shift error (rad)*/
    double frequency_m;          /**< Read in frequency of time varying field(Hz)*/

    bool fast_m;
    bool autophaseVeto_m;

    double designEnergy_m;

    Fieldmap* fieldmap_m;
    double startField_m;         /**< starting point of field(m)*/

private:
    double endField_m;

    CavityType type_m;

    static const boost::bimap<CavityType, std::string> bmCavityTypeString_s;

    double rmin_m;
    double rmax_m;
    double angle_m;
    double sinAngle_m;
    double cosAngle_m;
    double pdis_m;
    double gapwidth_m;
    double phi0_m;

    std::unique_ptr<double[]> RNormal_m;
    std::unique_ptr<double[]> VrNormal_m;
    std::unique_ptr<double[]> DvDr_m;
    int num_points_m;

    double getdE(const int& i,
                 const std::vector<double>& t,
                 const double& dz,
                 const double& phi,
                 const double& frequency,
                 const std::vector<double>& F) const;

    double getdT(const int& i,
                 const std::vector<double>& E,
                 const double& dz,
                 const double mass) const;

    double getdA(const int& i,
                 const std::vector<double>& t,
                 const double& dz,
                 const double& frequency,
                 const std::vector<double>& F) const;

    double getdB(const int& i,
                 const std::vector<double>& t,
                 const double& dz,
                 const double& frequency,
                 const std::vector<double>& F) const;

    // Not implemented.
    void operator=(const RFCavity&);
};

inline
double RFCavity::getdE(const int& i,
                       const std::vector<double>& t,
                       const double& dz,
                       const double& phi,
                       const double& frequency,
                       const std::vector<double>& F) const {
    return dz / (frequency * frequency * (t[i] - t[i-1]) * (t[i] - t[i-1])) *
        (frequency * (t[i] - t[i-1]) * (F[i] * std::sin(frequency * t[i] + phi) - F[i-1] * std::sin(frequency * t[i-1] + phi)) +
         (F[i] - F[i-1]) * (std::cos(frequency * t[i] + phi) - std::cos(frequency * t[i-1] + phi)));
}

inline
double RFCavity::getdT(const int& i,
                       const std::vector<double>& E,
                       const double& dz,
                       const double mass) const {
    double gamma1  = 1. + (19. * E[i-1] + 1. * E[i]) / (20. * mass);
    double gamma2  = 1. + (17. * E[i-1] + 3. * E[i]) / (20. * mass);
    double gamma3  = 1. + (15. * E[i-1] + 5. * E[i]) / (20. * mass);
    double gamma4  = 1. + (13. * E[i-1] + 7. * E[i]) / (20. * mass);
    double gamma5  = 1. + (11. * E[i-1] + 9. * E[i]) / (20. * mass);
    double gamma6  = 1. + (9. * E[i-1] + 11. * E[i]) / (20. * mass);
    double gamma7  = 1. + (7. * E[i-1] + 13. * E[i]) / (20. * mass);
    double gamma8  = 1. + (5. * E[i-1] + 15. * E[i]) / (20. * mass);
    double gamma9  = 1. + (3. * E[i-1] + 17. * E[i]) / (20. * mass);
    double gamma10 = 1. + (1. * E[i-1] + 19. * E[i]) / (20. * mass);
    return dz *
        (1. / std::sqrt(1. - 1. / (gamma1 * gamma1)) +
         1. / std::sqrt(1. - 1. / (gamma2 * gamma2)) +
         1. / std::sqrt(1. - 1. / (gamma3 * gamma3)) +
         1. / std::sqrt(1. - 1. / (gamma4 * gamma4)) +
         1. / std::sqrt(1. - 1. / (gamma5 * gamma5)) +
         1. / std::sqrt(1. - 1. / (gamma6 * gamma6)) +
         1. / std::sqrt(1. - 1. / (gamma7 * gamma7)) +
         1. / std::sqrt(1. - 1. / (gamma8 * gamma8)) +
         1. / std::sqrt(1. - 1. / (gamma9 * gamma9)) +
         1. / std::sqrt(1. - 1. / (gamma10 * gamma10))) / (10. * Physics::c);
}

inline
double RFCavity::getdA(const int& i,
                       const std::vector<double>& t,
                       const double& dz,
                       const double& frequency,
                       const std::vector<double>& F) const {
    double dt = t[i] - t[i-1];
    return dz / (frequency * frequency * dt * dt) *
        (frequency * dt * (F[i] * std::cos(frequency * t[i]) - F[i-1] * std::cos(frequency * t[i-1])) -
         (F[i] - F[i-1]) * (std::sin(frequency * t[i]) - std::sin(frequency * t[i-1])));
}

inline
double RFCavity::getdB(const int& i,
                       const std::vector<double>& t,
                       const double& dz,
                       const double& frequency,
                       const std::vector<double>& F) const {
    double dt = t[i] - t[i-1];
    return dz / (frequency * frequency * dt * dt) *
        (frequency * dt * (F[i] * std::sin(frequency * t[i]) - F[i-1] * std::sin(frequency * t[i-1])) +
         (F[i] - F[i-1]) * (std::cos(frequency * t[i]) - std::cos(frequency * t[i-1])));
}

inline
void RFCavity::setDesignEnergy(const double& ekin, bool) {
    designEnergy_m = ekin;
}

inline
double RFCavity::getDesignEnergy() const {
    return designEnergy_m;
}

inline
void RFCavity::dropFieldmaps() {
    fieldmap_m = nullptr;
}

inline
void RFCavity::setFieldMapFN(const std::string& fn) {
    filename_m = fn;
}

inline
void RFCavity::setAmplitudem(double vPeak) {
    scale_m = vPeak;
}

inline
double RFCavity::getAmplitudem() const {
    return scale_m;
}

inline
void RFCavity::setAmplitudeError(double vPeakError) {
    scaleError_m = vPeakError;
}

inline
double RFCavity::getAmplitudeError() const {
    return scaleError_m;
}

inline
void RFCavity::setFrequency(double freq) {
    frequency_m = freq;
}

inline
void RFCavity::setFrequencym(double freq) {
    frequency_m = freq;
}

inline
double RFCavity::getFrequencym() const {
    return frequency_m;
}

inline
void RFCavity::setPhasem(double phase) {
    phase_m = phase;
}

inline
double RFCavity::getPhasem() const {
    return phase_m;
}

inline
double RFCavity::getPhasem(double t) const {
    return phase_m + t * frequency_m;
}

inline
void RFCavity::setPhaseError(double phaseError) {
    phaseError_m = phaseError;
}

inline
double RFCavity::getPhaseError() const {
    return phaseError_m;
}

inline
CavityType RFCavity::getCavityType() const {
    return type_m;
}

inline
void RFCavity::setFast(bool fast) {
    fast_m = fast;
}

inline
bool RFCavity::getFast() const {
    return fast_m;
}

inline
void RFCavity::setAutophaseVeto(bool veto) {
    autophaseVeto_m = veto;
}

inline
bool RFCavity::getAutophaseVeto() const {
    return autophaseVeto_m;
}

inline
void RFCavity::setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> amplitudeTD) {
  amplitudeTD_m = amplitudeTD;
}

inline
void RFCavity::setAmplitudeModelName(std::string name) {
    amplitudeName_m=name;
}

inline
std::string RFCavity::getAmplitudeModelName() {
    return amplitudeName_m;
}

inline
void RFCavity::setPhaseModel(std::shared_ptr<AbstractTimeDependence> phaseTD) {
  phaseTD_m = phaseTD;
}

inline
void RFCavity::setPhaseModelName(std::string name) {
    phaseName_m=name;
}

inline
std::string RFCavity::getPhaseModelName() {
    return phaseName_m;
}

inline
void RFCavity::setFrequencyModel(std::shared_ptr<AbstractTimeDependence> frequencyTD) {
  frequencyTD_m = frequencyTD;
}

inline
void RFCavity::setFrequencyModelName(std::string name) {
  frequencyName_m=name;
}

inline
std::string RFCavity::getFrequencyModelName() {
    return frequencyName_m;
}

inline
CoordinateSystemTrafo RFCavity::getEdgeToBegin() const {
    CoordinateSystemTrafo ret(Vector_t(0, 0, startField_m),
                              Quaternion(1, 0, 0, 0));
    return ret;
}

inline
CoordinateSystemTrafo RFCavity::getEdgeToEnd() const {
    CoordinateSystemTrafo ret(Vector_t(0, 0, startField_m + getElementLength()),
                              Quaternion(1, 0, 0, 0));
    return ret;
}

#endif // CLASSIC_RFCavity_HH
