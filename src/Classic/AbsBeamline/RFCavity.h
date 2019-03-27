#ifndef CLASSIC_RFCavity_HH
#define CLASSIC_RFCavity_HH

// ------------------------------------------------------------------------
// $RCSfile: RFCavity.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFCavity
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------


#include "AbsBeamline/Component.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Physics/Physics.h"

class Fieldmap;

// Class RFCavity
// ------------------------------------------------------------------------
/// Interface for RF cavity.
//  Class RFCavity defines the abstract interface for RF cavities.


class RFCavity: public Component {

public:

    enum CavityType { SW, SGSW };
    /// Constructor with given name.
    explicit RFCavity(const std::string &name);

    RFCavity();
    RFCavity(const RFCavity &);
    virtual ~RFCavity();

    /// Apply visitor to RFCavity.
    virtual void accept(BeamlineVisitor &) const override;

    /// Get RF amplitude.
    virtual double getAmplitude() const = 0;

    /// Get RF frequencey.
    virtual double getFrequency() const = 0;

    /// Get RF phase.
    virtual double getPhase() const = 0;

    void dropFieldmaps();

    /// Set the name of the field map
    virtual void setFieldMapFN(std::string fmapfn);

    virtual std::string getFieldMapFN() const;

    virtual void setAmplitudem(double vPeak);
    virtual double getAmplitudem() const;
    virtual void setAmplitudeError(double vPeakError);
    virtual double getAmplitudeError() const;

    virtual void setFrequencym(double freq);

    void setFrequency(double freq);

    virtual double getFrequencym() const ;

    virtual void setPhasem(double phase);

    virtual double getPhasem() const;
    double getPhasem(double t) const;

    virtual void setPhaseError(double phaseError);
    virtual double getPhaseError() const;

    void setCavityType(std::string type);

    std::string getCavityType() const;

    virtual void setFast(bool fast);

    virtual bool getFast() const;

    virtual void setAutophaseVeto(bool veto = true);

    virtual bool getAutophaseVeto() const;

    virtual double getAutoPhaseEstimate(const double & E0, const double & t0, const double & q, const double & m);
    virtual double getAutoPhaseEstimateFallback(double E0, double t0, double q, double m);

    virtual std::pair<double, double> trackOnAxisParticle(const double & p0,
                                                          const double & t0,
                                                          const double & dt,
                                                          const double & q,
                                                          const double & mass,
							  std::ofstream *out = NULL);

    virtual void addKR(int i, double t, Vector_t &K) override;

    virtual void addKT(int i, double t, Vector_t &K) override;

    virtual bool apply(const size_t &i,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B) override;

    virtual bool apply(const Vector_t &R,
                       const Vector_t &P,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B) override;

    virtual bool applyToReferenceParticle(const Vector_t &R,
                                          const Vector_t &P,
                                          const double &t,
                                          Vector_t &E,
                                          Vector_t &B) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch,
                            std::shared_ptr<AbstractTimeDependence> freq_atd,
                            std::shared_ptr<AbstractTimeDependence> ampl_atd,
                            std::shared_ptr<AbstractTimeDependence> phase_atd);

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double &kineticEnergy) override;

    virtual void goOffline() override;

    virtual void setDesignEnergy(const double& ekin, bool changeable = true)  override;
    virtual double getDesignEnergy() const override;

    void setRmin(double rmin);

    void setRmax(double rmax);

    void setAzimuth(double angle);

    void setPerpenDistance(double pdis);

    void setGapWidth(double gapwidth);

    void setPhi0(double phi0);

    virtual double getRmin() const;

    virtual double getRmax() const;

    virtual double getAzimuth() const;

    virtual double getCosAzimuth() const;

    virtual double getSinAzimuth() const;

    virtual double getPerpenDistance() const;

    virtual double getGapWidth() const;

    virtual double getPhi0() const;

    virtual void setComponentType(std::string name) override;

    virtual std::string getComponentType()const override;

    virtual double getCycFrequency()const;

    void getMomentaKick(const double normalRadius, double momentum[], const double t, const double dtCorrt, const int PID, const double restMass,const int chargenumber);

    double spline(double z, double *za);

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    virtual bool isInside(const Vector_t &r) const override;

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
    virtual void getElementDimensions(double &begin,
                                      double &end) const override;

    virtual CoordinateSystemTrafo getEdgeToBegin() const override;
    virtual CoordinateSystemTrafo getEdgeToEnd() const override;

protected:
    std::shared_ptr<AbstractTimeDependence> phase_td_m;
    std::string phase_name_m;
    std::shared_ptr<AbstractTimeDependence> amplitude_td_m;
    std::string amplitude_name_m;
    std::shared_ptr<AbstractTimeDependence> frequency_td_m;
    std::string frequency_name_m;

    std::string filename_m;             /**< The name of the inputfile*/

    double scale_m;              /**< scale multiplier*/
    double scaleError_m;         /**< additive scale error*/
    double phase_m;              /**< phase shift of time varying field (rad)*/
    double phaseError_m;         /**< phase shift error (rad)*/
    double frequency_m;          /**< Read in frequency of time varying field(Hz)*/

    bool fast_m;
    bool autophaseVeto_m;

    double designEnergy_m;
private:
    Fieldmap* fieldmap_m;
    double startField_m;         /**< starting point of field(m)*/
    double endField_m;
    double length_m;

    CavityType type_m;

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

    double getdE(const int & i,
                 const std::vector<double> & t,
                 const double & dz,
                 const double & phi,
                 const double & frequency,
                 const std::vector<double> & F) const;

    double getdT(const int & i,
                 const std::vector<double> & E,
                 const double & dz,
                 const double mass) const;

    double getdA(const int & i,
                 const std::vector<double> & t,
                 const double & dz,
                 const double & frequency,
                 const std::vector<double> & F) const;

    double getdB(const int & i,
                 const std::vector<double> & t,
                 const double & dz,
                 const double & frequency,
                 const std::vector<double> & F) const;

    // Not implemented.
    void operator=(const RFCavity &);
};

inline
double RFCavity::getdE(const int & i,
                       const std::vector<double> & t,
                       const double & dz,
                       const double & phi,
                       const double & frequency,
                       const std::vector<double> & F) const {
    return dz / (frequency * frequency * (t[i] - t[i-1]) * (t[i] - t[i-1])) *
        (frequency * (t[i] - t[i-1]) * (F[i] * sin(frequency * t[i] + phi) - F[i-1] * sin(frequency * t[i-1] + phi)) +
         (F[i] - F[i-1]) * (cos(frequency * t[i] + phi) - cos(frequency * t[i-1] + phi)));
}

inline
double RFCavity::getdT(const int & i,
                       const std::vector<double> & E,
                       const double & dz,
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
        (1. / sqrt(1. - 1. / (gamma1 * gamma1)) +
         1. / sqrt(1. - 1. / (gamma2 * gamma2)) +
         1. / sqrt(1. - 1. / (gamma3 * gamma3)) +
         1. / sqrt(1. - 1. / (gamma4 * gamma4)) +
         1. / sqrt(1. - 1. / (gamma5 * gamma5)) +
         1. / sqrt(1. - 1. / (gamma6 * gamma6)) +
         1. / sqrt(1. - 1. / (gamma7 * gamma7)) +
         1. / sqrt(1. - 1. / (gamma8 * gamma8)) +
         1. / sqrt(1. - 1. / (gamma9 * gamma9)) +
         1. / sqrt(1. - 1. / (gamma10 * gamma10))) / (10. * Physics::c);
}

inline
double RFCavity::getdA(const int & i,
                       const std::vector<double> & t,
                       const double & dz,
                       const double & frequency,
                       const std::vector<double> & F) const {
    double dt = t[i] - t[i-1];
    return dz / (frequency * frequency * dt * dt) *
        (frequency * dt * (F[i] * cos(frequency * t[i]) - F[i-1] * cos(frequency * t[i-1])) -
         (F[i] - F[i-1]) * (sin(frequency * t[i]) - sin(frequency * t[i-1])));
}

inline
double RFCavity::getdB(const int & i,
                       const std::vector<double> & t,
                       const double & dz,
                       const double & frequency,
                       const std::vector<double> & F) const {
    double dt = t[i] - t[i-1];
    return dz / (frequency * frequency * dt * dt) *
        (frequency * dt * (F[i] * sin(frequency * t[i]) - F[i-1] * sin(frequency * t[i-1])) +
         (F[i] - F[i-1]) * (cos(frequency * t[i]) - cos(frequency * t[i-1])));
}

inline
void RFCavity::setDesignEnergy(const double& ekin, bool)
{
    designEnergy_m = ekin;
}

inline
double RFCavity::getDesignEnergy() const
{
    return designEnergy_m;
}

inline
void RFCavity::dropFieldmaps() {
    fieldmap_m = NULL;
}

inline
void RFCavity::setFieldMapFN(std::string fn) {
    filename_m = fn;
}

inline
std::string RFCavity::getFieldMapFN() const {
    return filename_m;
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
void RFCavity::setCavityType(std::string type) {

}

inline
std::string RFCavity::getCavityType() const {
    return "SW";
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
void RFCavity::setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> amplitude_td) {
  amplitude_td_m = amplitude_td;
}

inline
void RFCavity::setAmplitudeModelName(std::string name) {
    amplitude_name_m=name;
}

inline
std::string RFCavity::getAmplitudeModelName() {
    return amplitude_name_m;
}

inline
void RFCavity::setPhaseModel(std::shared_ptr<AbstractTimeDependence> phase_td) {
  phase_td_m = phase_td;
}

inline
void RFCavity::setPhaseModelName(std::string name) {
    phase_name_m=name;
}

inline
std::string RFCavity::getPhaseModelName() {
    return phase_name_m;
}

inline
void RFCavity::setFrequencyModel(std::shared_ptr<AbstractTimeDependence> frequency_td) {
  frequency_td_m = frequency_td;
}

inline
void RFCavity::setFrequencyModelName(std::string name) {
  frequency_name_m=name;
}

inline
std::string RFCavity::getFrequencyModelName() {
    return frequency_name_m;
}


inline
CoordinateSystemTrafo RFCavity::getEdgeToBegin() const
{
    CoordinateSystemTrafo ret(Vector_t(0, 0, startField_m),
                              Quaternion(1, 0, 0, 0));

    return ret;
}

inline
CoordinateSystemTrafo RFCavity::getEdgeToEnd() const
{
    CoordinateSystemTrafo ret(Vector_t(0, 0, endField_m),
                              Quaternion(1, 0, 0, 0));

    return ret;
}


#endif // CLASSIC_RFCavity_HH
