//
// Class TravelingWave
//   Defines the abstract interface for Traveling Wave.
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
#ifndef CLASSIC_TravelingWave_HH
#define CLASSIC_TravelingWave_HH

#include "AbsBeamline/RFCavity.h"
#include "Physics/Physics.h"

#include <cmath>

class Fieldmap;

class TravelingWave: public RFCavity {

public:

    /// Constructor with given name.
    explicit TravelingWave(const std::string& name);

    TravelingWave();
    TravelingWave(const TravelingWave&);
    virtual ~TravelingWave();

    /// Apply visitor to TravelingWave.
    virtual void accept(BeamlineVisitor&) const override;

    /// Get RF amplitude.
    virtual double getAmplitude() const override = 0;

    /// Get RF frequency.
    virtual double getFrequency() const override = 0;

    /// Get RF phase.
    virtual double getPhase() const override = 0;

    virtual void setPhasem(double phase) override;

    void setNumCells(int NumCells);

    void setMode(double mode);

    virtual double getAutoPhaseEstimate(const double& E0, const double& t0, const double& q, const double& m) override;

    virtual bool apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B) override;

    virtual bool apply(const Vector_t& R, const Vector_t& P, const double& t, Vector_t& E, Vector_t& B) override;

    virtual bool applyToReferenceParticle(const Vector_t& R, const Vector_t& P, const double& t, Vector_t& E, Vector_t& B) override;

    virtual void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double& kineticEnergy) override;

    virtual void goOffline() override;

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual bool isInside(const Vector_t& r) const override;

    virtual void getElementDimensions(double& begin,
                                      double& end) const override;

    virtual CoordinateSystemTrafo getEdgeToBegin() const override;
    virtual CoordinateSystemTrafo getEdgeToEnd() const override;

private:
    double scaleCore_m;
    double scaleCoreError_m;

    double phaseCore1_m;
    double phaseCore2_m;
    double phaseExit_m;

    double startCoreField_m;         /**< starting point of field(m)*/
    double startExitField_m;
    double mappedStartExitField_m;

    double periodLength_m;
    int numCells_m;
    double cellLength_m;
    double mode_m;

    inline double getdE(const int& i,
                        const int& I,
                        const std::vector<double>& t,
                        const double& phi,
                        const std::vector<std::pair<double, double> >& F) const;

    inline double getdT(const int& i,
                        const int& I,
                        const std::vector<double>& E,
                        const std::vector<std::pair<double, double> >& F,
                        const double mass) const;

    inline double getdA(const int& i,
                        const int& I,
                        const std::vector<double>& t,
                        const double& phi,
                        const std::vector<std::pair<double, double> >& F) const;

    inline double getdB(const int& i,
                        const int& I,
                        const std::vector<double>& t,
                        const double& phi,
                        const std::vector<std::pair<double, double> >& F) const;
    // Not implemented.
    void operator=(const TravelingWave&);
};

double TravelingWave::getdE(const int& i,
                            const int& I,
                            const std::vector<double>& t,
                            const double& phi,
                            const std::vector<std::pair<double, double> >& F) const {
    return (F[I].first - F[I-1].first) / (frequency_m * frequency_m * (t[i] - t[i-1]) * (t[i] - t[i-1])) *
        (frequency_m * (t[i] - t[i-1]) * (F[I].second * std::sin(frequency_m * t[i] + phi) - F[I-1].second * std::sin(frequency_m * t[i-1] + phi)) +
         (F[I].second - F[I-1].second) * (std::cos(frequency_m * t[i] + phi) - std::cos(frequency_m * t[i-1] + phi)));
}

double TravelingWave::getdT(const int& i,
                            const int& I,
                            const std::vector<double>& E,
                            const std::vector<std::pair<double, double> >& F,
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
    return (F[I].first - F[I-1].first) *
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

double TravelingWave::getdA(const int& i,
                            const int& I,
                            const std::vector<double>& t,
                            const double& phi,
                            const std::vector<std::pair<double, double> >& F) const {
    double dt = t[i] - t[i-1];
    return (F[I].first - F[I-1].first) / (frequency_m * frequency_m * dt * dt) *
        (frequency_m * dt * (F[I].second * std::cos(frequency_m * t[i] + phi) - F[I-1].second * std::cos(frequency_m * t[i-1] + phi)) -
         (F[I].second - F[I-1].second) * (std::sin(frequency_m * t[i] + phi) - std::sin(frequency_m * t[i-1] + phi)));
}

double TravelingWave::getdB(const int& i,
                            const int& I,
                            const std::vector<double>& t,
                            const double& phi,
                            const std::vector<std::pair<double, double> >& F) const {
    double dt = t[i] - t[i-1];
    return (F[I].first - F[I-1].first) / (frequency_m * frequency_m * dt * dt) *
        (frequency_m * dt * (F[I].second * std::sin(frequency_m * t[i] + phi) - F[I-1].second * std::sin(frequency_m * t[i-1] + phi)) +
         (F[I].second - F[I-1].second) * (std::cos(frequency_m * t[i] + phi) - std::cos(frequency_m * t[i-1] + phi)));
}

inline
void TravelingWave::setPhasem(double phase) {
    phase_m = phase;
    phaseCore1_m = phase_m + Physics::pi * mode_m / 2.0;
    phaseCore2_m = phase_m + Physics::pi * mode_m * 1.5;
    phaseExit_m = phase_m - Physics::two_pi * ((numCells_m - 1) * mode_m - std::floor((numCells_m - 1) * mode_m));
}

inline
void TravelingWave::setNumCells(int NumCells) {
    numCells_m = NumCells;
}

inline
void TravelingWave::setMode(double mode) {
    mode_m = mode;
}

inline
CoordinateSystemTrafo TravelingWave::getEdgeToBegin() const {
    CoordinateSystemTrafo ret(Vector_t(0, 0, -0.5 * periodLength_m),
                              Quaternion(1, 0, 0, 0));
    return ret;
}

inline
CoordinateSystemTrafo TravelingWave::getEdgeToEnd() const {
    CoordinateSystemTrafo ret(Vector_t(0, 0, -0.5 * periodLength_m + getElementLength()),
                              Quaternion(1, 0, 0, 0));
    return ret;
}

#endif // CLASSIC_TravelingWave_HH
