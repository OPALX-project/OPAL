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
#include "AbsBeamline/RFCavity.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Units.h"
#include "Steppers/BorisPusher.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

const boost::bimap<CavityType, std::string> RFCavity::bmCavityTypeString_s =
    boost::assign::list_of<const boost::bimap<CavityType, std::string>::relation>
        (CavityType::SW,   "STANDING")
        (CavityType::SGSW, "SINGLEGAP");


RFCavity::RFCavity():
    RFCavity("")
{}

RFCavity::RFCavity(const RFCavity& right):
    Component(right),
    phaseTD_m(right.phaseTD_m),
    phaseName_m(right.phaseName_m),
    amplitudeTD_m(right.amplitudeTD_m),
    amplitudeName_m(right.amplitudeName_m),
    frequencyTD_m(right.frequencyTD_m),
    frequencyName_m(right.frequencyName_m),
    filename_m(right.filename_m),
    scale_m(right.scale_m),
    scaleError_m(right.scaleError_m),
    phase_m(right.phase_m),
    phaseError_m(right.phaseError_m),
    frequency_m(right.frequency_m),
    fast_m(right.fast_m),
    autophaseVeto_m(right.autophaseVeto_m),
    designEnergy_m(right.designEnergy_m),
    fieldmap_m(right.fieldmap_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    type_m(right.type_m),
    rmin_m(right.rmin_m),
    rmax_m(right.rmax_m),
    angle_m(right.angle_m),
    sinAngle_m(right.sinAngle_m),
    cosAngle_m(right.cosAngle_m),
    pdis_m(right.pdis_m),
    gapwidth_m(right.gapwidth_m),
    phi0_m(right.phi0_m),
    RNormal_m(nullptr),
    VrNormal_m(nullptr),
    DvDr_m(nullptr),
    num_points_m(right.num_points_m)
{}

RFCavity::RFCavity(const std::string& name):
    Component(name),
    phaseTD_m(nullptr),
    amplitudeTD_m(nullptr),
    frequencyTD_m(nullptr),
    filename_m(""),
    scale_m(1.0),
    scaleError_m(0.0),
    phase_m(0.0),
    phaseError_m(0.0),
    frequency_m(0.0),
    fast_m(true),
    autophaseVeto_m(false),
    designEnergy_m(-1.0),
    fieldmap_m(nullptr),
    startField_m(0.0),
    endField_m(0.0),
    type_m(CavityType::SW),
    rmin_m(0.0),
    rmax_m(0.0),
    angle_m(0.0),
    sinAngle_m(0.0),
    cosAngle_m(0.0),
    pdis_m(0.0),
    gapwidth_m(0.0),
    phi0_m(0.0),
    RNormal_m(nullptr),
    VrNormal_m(nullptr),
    DvDr_m(nullptr),
    num_points_m(0)
{}

RFCavity::~RFCavity() {
}

void RFCavity::accept(BeamlineVisitor& visitor) const {
    visitor.visitRFCavity(*this);
}

bool RFCavity::apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool RFCavity::apply(const Vector_t& R,
                     const Vector_t& /*P*/,
                     const double& t,
                     Vector_t& E,
                     Vector_t& B) {

    if (R(2) >= startField_m &&
        R(2) < startField_m + getElementLength()) {
        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        bool outOfBounds = fieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) return getFlagDeleteOnTransverseExit();

        E += (scale_m + scaleError_m) * std::cos(frequency_m * t + phase_m + phaseError_m) * tmpE;
        B -= (scale_m + scaleError_m) * std::sin(frequency_m * t + phase_m + phaseError_m) * tmpB;

    }
    return false;
}

bool RFCavity::applyToReferenceParticle(const Vector_t& R,
                                        const Vector_t& /*P*/,
                                        const double& t,
                                        Vector_t& E,
                                        Vector_t& B) {

    if (R(2) >= startField_m &&
        R(2) < startField_m + getElementLength()) {
        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        bool outOfBounds = fieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) return true;

        E += scale_m * std::cos(frequency_m * t + phase_m) * tmpE;
        B -= scale_m * std::sin(frequency_m * t + phase_m) * tmpB;

    }
    return false;
}

void RFCavity::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {

    startField_m = endField_m = 0.0;
    if (bunch == nullptr) {
        startField = startField_m;
        endField = endField_m;

        return;
    }

    Inform msg("RFCavity ", *gmsg);
    std::stringstream errormsg;
    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    fieldmap_m->getFieldDimensions(startField_m, endField);
    if (endField <= startField_m) {
        throw GeneralClassicException("RFCavity::initialise",
                                      "The length of the field map '" + filename_m + "' is zero or negative");
    }

    msg << level2 << getName() << " using file ";
    fieldmap_m->getInfo(&msg);
    if (std::abs((frequency_m - fieldmap_m->getFrequency()) / frequency_m) > 0.01) {
        errormsg << "FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP '" << filename_m << "';\n"
                 << frequency_m / Physics::two_pi * Units::Hz2MHz << " MHz <> "
                 << fieldmap_m->getFrequency() / Physics::two_pi * Units::Hz2MHz << " MHz; TAKE ON THE LATTER";
        std::string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
        ERRORMSG(errormsg_str << "\n" << endl);
        if (Ippl::myNode() == 0) {
            std::ofstream omsg("errormsg.txt", std::ios_base::app);
            omsg << errormsg_str << std::endl;
            omsg.close();
        }
        frequency_m = fieldmap_m->getFrequency();
    }
    setElementLength(endField - startField_m);
}

// In current version ,this function reads in the cavity voltage profile data from file.
void RFCavity::initialise(PartBunchBase<double, 3>* bunch,
                          std::shared_ptr<AbstractTimeDependence> freq_atd,
                          std::shared_ptr<AbstractTimeDependence> ampl_atd,
                          std::shared_ptr<AbstractTimeDependence> phase_atd) {

    RefPartBunch_m = bunch;

    /// set the time dependent models
    setAmplitudeModel(ampl_atd);
    setPhaseModel(phase_atd);
    setFrequencyModel(freq_atd);

    std::ifstream in(filename_m.c_str());
    in >> num_points_m;

    RNormal_m  = std::unique_ptr<double[]>(new double[num_points_m]);
    VrNormal_m = std::unique_ptr<double[]>(new double[num_points_m]);
    DvDr_m     = std::unique_ptr<double[]>(new double[num_points_m]);

    for (int i = 0; i < num_points_m; i++) {
        if (in.eof()) {
            throw GeneralClassicException("RFCavity::initialise",
                                          "Not enough data in file '" + filename_m +
                                          "', please check the data format");
        }
        in >> RNormal_m[i] >> VrNormal_m[i] >> DvDr_m[i];

        VrNormal_m[i] *= RefPartBunch_m->getQ();
        DvDr_m[i]     *= RefPartBunch_m->getQ();
    }
    sinAngle_m = std::sin(angle_m * Units::deg2rad);
    cosAngle_m = std::cos(angle_m * Units::deg2rad);

    if (!frequencyName_m.empty()) {
      *gmsg << "* Timedependent frequency model " << frequencyName_m << endl;
    }

    *gmsg << "* Cavity voltage data read successfully!" << endl;
}

void RFCavity::finalise()
{}

bool RFCavity::bends() const {
    return false;
}

void RFCavity::goOnline(const double&) {
    Fieldmap::readMap(filename_m);

    online_m = true;
}

void RFCavity::goOffline() {
    Fieldmap::freeMap(filename_m);

    online_m = false;
}

void  RFCavity::setRmin(double rmin) {
    rmin_m = rmin;
}

void  RFCavity::setRmax(double rmax) {
    rmax_m = rmax;
}

void  RFCavity::setAzimuth(double angle) {
    angle_m = angle;
}

void  RFCavity::setPerpenDistance(double pdis) {
    pdis_m = pdis;
}

void  RFCavity::setGapWidth(double gapwidth) {
    gapwidth_m = gapwidth;
}

void RFCavity::setPhi0(double phi0) {
    phi0_m = phi0;
}

double  RFCavity::getRmin() const {
    return rmin_m;
}

double  RFCavity::getRmax() const {
    return rmax_m;
}

double  RFCavity::getAzimuth() const {
    return angle_m;
}

double  RFCavity::getSinAzimuth() const {
    return sinAngle_m;
}

double  RFCavity::getCosAzimuth() const {
    return cosAngle_m;
}

double  RFCavity::getPerpenDistance() const {
    return pdis_m;
}

double  RFCavity::getGapWidth() const {
    return gapwidth_m;
}

double RFCavity::getPhi0() const {
    return phi0_m;
}

void RFCavity::setCavityType(const std::string& name) {
    auto it = bmCavityTypeString_s.right.find(name);
    if (it != bmCavityTypeString_s.right.end()) {
        type_m = it->second;
    } else {
        type_m = CavityType::SW;
    }
}

std::string RFCavity::getCavityTypeString() const {
    return bmCavityTypeString_s.left.at(type_m);
}

std::string RFCavity::getFieldMapFN() const {
    if (filename_m.empty()) {
        throw GeneralClassicException("RFCavity::getFieldMapFN",
                                      "The attribute \"FMAPFN\" isn't set "
                                      "for the \"RFCAVITY\" element!");
    } else if (boost::filesystem::exists(filename_m)) {
        return filename_m;
    } else {
        throw GeneralClassicException("RFCavity::getFieldMapFN",
                                      "Failed to open file '" + filename_m +
                                      "', please check if it exists");
    }
}

double RFCavity::getCycFrequency() const {
    return frequency_m;
}

/**
   \brief used in OPAL-cycl

   Is called from OPAL-cycl and can handle
   time dependent frequency, amplitude and phase

   At the moment (test) only the frequency is time
   dependent

 */
void RFCavity::getMomentaKick(const double normalRadius,
                              double momentum[],
                              const double t,
                              const double dtCorrt,
                              const int PID,
                              const double restMass,
                              const int chargenumber) {

    double derivate;

    double momentum2  = momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2];
    double betgam = std::sqrt(momentum2);

    double gamma = std::sqrt(1.0 + momentum2);
    double beta = betgam / gamma;

    double Voltage = spline(normalRadius, &derivate) * scale_m * Units::MVpm2Vpm;

    double Ufactor = 1.0;

    double frequency = frequency_m * frequencyTD_m->getValue(t);

    if (gapwidth_m > 0.0) {
    	double transit_factor = 0.5 * frequency * gapwidth_m * Units::mm2m / (Physics::c * beta);
        Ufactor = std::sin(transit_factor) / transit_factor;
    }

    Voltage *= Ufactor;
    // rad/s, ns --> rad
    double nphase = (frequency * (t + dtCorrt) * Units::ns2s) - phi0_m * Units::deg2rad;
    double dgam = Voltage * std::cos(nphase) / (restMass);

    double tempdegree = std::fmod(nphase * Units::rad2deg, 360.0);
    if (tempdegree > 270.0) tempdegree -= 360.0;

    gamma += dgam;

    double newmomentum2 = std::pow(gamma, 2) - 1.0;

    double pr = momentum[0] * cosAngle_m + momentum[1] * sinAngle_m;
    double ptheta = std::sqrt(newmomentum2 - std::pow(pr, 2));
    double px = pr * cosAngle_m - ptheta * sinAngle_m ; // x
    double py = pr * sinAngle_m + ptheta * cosAngle_m; // y

    double rotate = -derivate * (scale_m * Units::MVpm2Vpm) / ((rmax_m - rmin_m) * Units::mm2m) * std::sin(nphase) / (frequency * Physics::two_pi) / (betgam * restMass / Physics::c / chargenumber); // radian

    /// B field effects
    momentum[0] =  std::cos(rotate) * px + std::sin(rotate) * py;
    momentum[1] = -std::sin(rotate) * px + std::cos(rotate) * py;

    if (PID == 0) {
        Inform  m("OPAL", *gmsg, Ippl::myNode());

        m << "* Cavity " << getName() << " Phase= " << tempdegree << " [deg] transit time factor=  " << Ufactor
          << " dE= " << dgam *restMass * Units::eV2MeV << " [MeV]"
          << " E_kin= " << (gamma - 1.0)*restMass * Units::eV2MeV << " [MeV] Time dep freq = " << frequencyTD_m->getValue(t) << endl;
    }

}

/* cubic spline subrutine */
double RFCavity::spline(double z, double *za) {
    double splint;

    // domain-test and handling of case "1-support-point"
    if (num_points_m < 1) {
        throw GeneralClassicException("RFCavity::spline",
                                      "no support points!");
    }
    if (num_points_m == 1) {
        splint = RNormal_m[0];
        *za = 0.0;
        return splint;
    }

    // search the two support-points
    int il, ih;
    il = 0;
    ih = num_points_m - 1;
    while ((ih - il) > 1) {
        int i = (int)((il + ih) / 2.0);
        if (z < RNormal_m[i]) {
            ih = i;
        } else if (z > RNormal_m[i]) {
            il = i;
        } else if (z == RNormal_m[i]) {
            il = i;
            ih = i + 1;
            break;
        }
    }

    double x1 =  RNormal_m[il];
    double x2 =  RNormal_m[ih];
    double y1 =  VrNormal_m[il];
    double y2 =  VrNormal_m[ih];
    double y1a = DvDr_m[il];
    double y2a = DvDr_m[ih];
    //
    // determination of the requested function-values and its derivatives
    //
    double dx  = x2 - x1;
    double dy  = y2 - y1;
    double u   = (z - x1) / dx;
    double u2  = u * u;
    double u3  = u2 * u;
    double dy2 = -2.0 * dy;
    double ya2 = y2a + 2.0 * y1a;
    double dy3 = 3.0 * dy;
    double ya3 = y2a + y1a;
    double yb2 = dy2 + dx * ya3;
    double yb4 = dy3 - dx * ya2;
    splint = y1  + u * dx * y1a +       u2 * yb4 +        u3 * yb2;
    *za    =            y1a + 2.0 * u / dx * yb4 + 3.0 * u2 / dx * yb2;
    // if(m>=1) za=y1a+2.0*u/dx*yb4+3.0*u2/dx*yb2;
    // if(m>=2) za[1]=2.0/dx2*yb4+6.0*u/dx2*yb2;
    // if(m>=3) za[2]=6.0/dx3*yb2;

    return splint;
}

void RFCavity::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


ElementType RFCavity::getType() const {
    return ElementType::RFCAVITY;
}

double RFCavity::getAutoPhaseEstimateFallback(double E0, double t0, double q, double mass) {

    const double dt = 1e-13;
    const double p0 = Util::getBetaGamma(E0, mass);
    const double origPhase =getPhasem();
    double dphi = Physics::pi / 18;

    double phi = 0.0;
    setPhasem(phi);
    std::pair<double, double> ret = trackOnAxisParticle(p0, t0, dt, q, mass);
    double phimax = 0.0;
    double Emax = Util::getKineticEnergy(Vector_t(0.0, 0.0, ret.first), mass);
    phi += dphi;

    for (unsigned int j = 0; j < 2; ++ j) {
        for (unsigned int i = 0; i < 36; ++ i, phi += dphi) {
            setPhasem(phi);
            ret = trackOnAxisParticle(p0, t0, dt, q, mass);
            double Ekin = Util::getKineticEnergy(Vector_t(0.0, 0.0, ret.first), mass);
            if (Ekin > Emax) {
                Emax = Ekin;
                phimax = phi;
            }
        }

        phi = phimax - dphi;
        dphi = dphi / 17.5;
    }

    phimax = phimax - std::round(phimax / Physics::two_pi) * Physics::two_pi;
    phimax = std::fmod(phimax, Physics::two_pi);

    const int prevPrecision = Ippl::Info->precision(8);
    INFOMSG(level2
            << "estimated phase= " << phimax << " rad = "
            << phimax * Units::rad2deg << " deg \n"
            << "Ekin= " << Emax << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);

    setPhasem(origPhase);
    return phimax;
}

double RFCavity::getAutoPhaseEstimate(const double& E0, const double& t0,
                                      const double& q, const double& mass) {
    std::vector<double> t, E, t2, E2;
    std::vector<double> F;
    std::vector< std::pair< double, double > > G;
    gsl_spline *onAxisInterpolants;
    gsl_interp_accel *onAxisAccel;

    double phi = 0.0, tmp_phi, dphi = 0.5 * Units::deg2rad;
    double dz = 1.0, length = 0.0;
    fieldmap_m->getOnaxisEz(G);
    if (G.size() == 0) return 0.0;
    double begin = (G.front()).first;
    double end   = (G.back()).first;
    std::unique_ptr<double[]> zvals(      new double[G.size()]);
    std::unique_ptr<double[]> onAxisField(new double[G.size()]);

    for (size_t j = 0; j < G.size(); ++ j) {
        zvals[j] = G[j].first;
        onAxisField[j] = G[j].second;
    }
    onAxisInterpolants = gsl_spline_alloc(gsl_interp_cspline, G.size());
    onAxisAccel = gsl_interp_accel_alloc();
    gsl_spline_init(onAxisInterpolants, zvals.get(), onAxisField.get(), G.size());

    length = end - begin;
    dz = length / G.size();

    G.clear();

    unsigned int N = (int)std::floor(length / dz + 1);
    dz = length / N;

    F.resize(N);
    double z = begin;
    for (size_t j = 0; j < N; ++ j, z += dz) {
        F[j] = gsl_spline_eval(onAxisInterpolants, z, onAxisAccel);
    }
    gsl_spline_free(onAxisInterpolants);
    gsl_interp_accel_free(onAxisAccel);

    t.resize(N, t0);
    t2.resize(N, t0);
    E.resize(N, E0);
    E2.resize(N, E0);

    z = begin + dz;
    for (unsigned int i = 1; i < N; ++ i, z += dz) {
        E[i] = E[i - 1] + dz * scale_m / mass;
        E2[i] = E[i];
    }

    for (int iter = 0; iter < 10; ++ iter) {
        double A = 0.0;
        double B = 0.0;
        for (unsigned int i = 1; i < N; ++ i) {
            t[i] = t[i - 1] + getdT(i, E, dz, mass);
            t2[i] = t2[i - 1] + getdT(i, E2, dz, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, t, dz, frequency_m, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, t, dz, frequency_m, F);
        }

        if (std::abs(B) > 0.0000001) {
            tmp_phi = std::atan(A / B);
        } else {
            tmp_phi = Physics::pi / 2;
        }
        if (q * (A * std::sin(tmp_phi) + B * std::cos(tmp_phi)) < 0) {
            tmp_phi += Physics::pi;
        }

        if (std::abs (phi - tmp_phi) < frequency_m * (t[N - 1] - t[0]) / (10 * N)) {
            for (unsigned int i = 1; i < N; ++ i) {
                E[i] = E[i - 1];
                E[i] += q * scale_m * getdE(i, t, dz, phi, frequency_m, F) ;
            }
            const int prevPrecision = Ippl::Info->precision(8);
            INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad = "
                    << tmp_phi * Units::rad2deg << " deg \n"
                    << "Ekin= " << E[N - 1] << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);

            return tmp_phi;
        }
        phi = tmp_phi - std::round(tmp_phi / Physics::two_pi) * Physics::two_pi;

        for (unsigned int i = 1; i < N; ++ i) {
            E[i] = E[i - 1];
            E2[i] = E2[i - 1];
            E[i] += q * scale_m * getdE(i, t, dz, phi, frequency_m, F) ;
            E2[i] += q * scale_m * getdE(i, t2, dz, phi + dphi, frequency_m, F);
            double a = E[i], b = E2[i];
            if (std::isnan(a) || std::isnan(b)) {
                return getAutoPhaseEstimateFallback(E0, t0, q, mass);
            }
            t[i] = t[i - 1] + getdT(i, E, dz, mass);
            t2[i] = t2[i - 1] + getdT(i, E2, dz, mass);

            E[i]  = E [i - 1];
            E2[i] = E2[i - 1];
            E[i]  += q * scale_m * getdE(i, t, dz, phi, frequency_m, F) ;
            E2[i] += q * scale_m * getdE(i, t2, dz, phi + dphi, frequency_m, F);
        }

        double cosine_part = 0.0, sine_part = 0.0;
        double p0 = Util::getBetaGamma(E0, mass);
        cosine_part += scale_m * std::cos(frequency_m * t0) * F[0];
        sine_part += scale_m * std::sin(frequency_m * t0) * F[0];

        double totalEz0 = std::cos(phi) * cosine_part - std::sin(phi) * sine_part;

        if (p0 + q * totalEz0 * (t[1] - t[0]) * Physics::c / mass < 0) {
            // make totalEz0 = 0
            tmp_phi = std::atan(cosine_part / sine_part);
            if (std::abs (tmp_phi - phi) > Physics::pi) {
                phi = tmp_phi + Physics::pi;
            } else {
                phi = tmp_phi;
            }
        }
    }

    const int prevPrecision = Ippl::Info->precision(8);
    INFOMSG(level2
            << "estimated phase= " << tmp_phi << " rad = "
            << tmp_phi * Units::rad2deg << " deg \n"
            << "Ekin= " << E[N - 1] << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);

    return phi;
}

std::pair<double, double> RFCavity::trackOnAxisParticle(const double& p0,
                                                        const double& t0,
                                                        const double& dt,
                                                        const double& /*q*/,
                                                        const double& mass,
                                                        std::ofstream *out) {
    Vector_t p(0, 0, p0);
    double t = t0;
    BorisPusher integrator(*RefPartBunch_m->getReference());
    const double cdt = Physics::c * dt;
    const double zbegin = startField_m;
    const double zend = getElementLength() + startField_m;

    Vector_t z(0.0, 0.0, zbegin);
    double dz = 0.5 * p(2) / Util::getGamma(p) * cdt;
    Vector_t Ef(0.0), Bf(0.0);

    if (out) *out << std::setw(18) << z[2]
                  << std::setw(18) << Util::getKineticEnergy(p, mass)
                  << std::endl;
    while (z(2) + dz < zend && z(2) + dz > zbegin) {
        z /= cdt;
        integrator.push(z, p, dt);
        z *= cdt;

        Ef = 0.0;
        Bf = 0.0;
        if (z(2) >= zbegin && z(2) <= zend) {
            applyToReferenceParticle(z, p, t + 0.5 * dt, Ef, Bf);
        }
        integrator.kick(z, p, Ef, Bf, dt);

        dz = 0.5 * p(2) / std::sqrt(1.0 + dot(p, p)) * cdt;
        z /= cdt;
        integrator.push(z, p, dt);
        z *= cdt;
        t += dt;

        if (out) *out << std::setw(18) << z[2]
                      << std::setw(18) << Util::getKineticEnergy(p, mass)
                      << std::endl;
    }

    const double beta = std::sqrt(1. - 1 / (dot(p, p) + 1.));
    const double tErr = (z(2) - zend) / (Physics::c * beta);

    return std::pair<double, double>(p(2), t - tErr);
}

bool RFCavity::isInside(const Vector_t& r) const {
    if (isInsideTransverse(r)) {
        return fieldmap_m->isInside(r);
    }

    return false;
}

double RFCavity::getElementLength() const {
    double length = ElementBase::getElementLength();
    if (length < 1e-10 && fieldmap_m != nullptr) {
        double start, end;
        fieldmap_m->getFieldDimensions(start, end);
        length = end - start;
    }

    return length;
}

void RFCavity::getElementDimensions(double& begin, double& end) const {
    fieldmap_m->getFieldDimensions(begin, end);
}
