// ------------------------------------------------------------------------
// $RCSfile: RFCavity.cpp,v $
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
// $Author: andreas adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Steppers/BorisPusher.h"
#include "Fields/Fieldmap.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

using namespace std;

// Class RFCavity
// ------------------------------------------------------------------------

RFCavity::RFCavity():
    RFCavity("")
{}


RFCavity::RFCavity(const RFCavity &right):
    Component(right),
    phase_td_m(right.phase_td_m),
    phase_name_m(right.phase_name_m),
    amplitude_td_m(right.amplitude_td_m),
    amplitude_name_m(right.amplitude_name_m),
    frequency_td_m(right.frequency_td_m),
    frequency_name_m(right.frequency_name_m),
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
    length_m(right.length_m),
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
{
    setElType(isRF);
}


RFCavity::RFCavity(const std::string &name):
    Component(name),
    phase_td_m(nullptr),
    amplitude_td_m(nullptr),
    frequency_td_m(nullptr),
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
    length_m(0.0),
    type_m(SW),
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
    //     RNormal_m(std::nullptr_t(NULL)),
    //     VrNormal_m(std::nullptr_t(NULL)),
    //     DvDr_m(std::nullptr_t(NULL)),
    num_points_m(0)
{
    setElType(isRF);
}


RFCavity::~RFCavity() {
    // FIXME: in deleteFielmak, a map find makes problems
    //       Fieldmap::deleteFieldmap(filename_m);
    //~ if(RNormal_m) {
    //~ delete[] RNormal_m;
    //~ delete[] VrNormal_m;
    //~ delete[] DvDr_m;
    //~ }
}

void RFCavity::accept(BeamlineVisitor &visitor) const {
    visitor.visitRFCavity(*this);
}

/**
 * ENVELOPE COMPONENT for radial focussing of the beam
 * Calculates the transverse envelope component for the RF cavity
 * element and adds it to the K vector
*/
void RFCavity::addKR(int i, double t, Vector_t &K) {

    double pz = RefPartBunch_m->getZ(i) - startField_m;
    if (pz < 0.0 ||
        pz >= length_m) return;

    const Vector_t tmpR(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), pz);
    double k = -Physics::q_e / (2.0 * RefPartBunch_m->getGamma(i) * Physics::EMASS);

    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    double Ez = tmpE(2);

    tmpE = Vector_t(0.0);
    fieldmap_m->getFieldDerivative(tmpR, tmpE, tmpB, DZ);

    double wtf = frequency_m * t + phase_m;
    double kj = k * scale_m * (tmpE(2) * cos(wtf) - RefPartBunch_m->getBeta(i) * frequency_m * Ez * sin(wtf) / Physics::c);

    K += Vector_t(kj, kj, 0.0);
}


/**
 * ENVELOPE COMPONENT for transverse kick (only has an impact if x0, y0 != 0)
 * Calculates the transverse kick component for the RF cavity element and adds it to
 * the K vector. Only important for off track tracking, otherwise KT = 0.
*/
void RFCavity::addKT(int i, double t, Vector_t &K) {

    RefPartBunch_m->actT();

    double pz = RefPartBunch_m->getZ(i) - startField_m;
    if (pz < 0.0 ||
        pz >= length_m) return;

    //XXX: BET parameter, default is 1.
    //If cxy != 1, then cxy = true
    bool cxy = false; // default
    double kx = 0.0, ky = 0.0;
    if(cxy) {
        const Vector_t tmpA(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), pz);

        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
        fieldmap_m->getFieldstrength(tmpA, tmpE, tmpB);

        double cwtf = cos(frequency_m * t + phase_m);
        double cf = -Physics::q_e / (RefPartBunch_m->getGamma(i) * Physics::m_e);
        kx += -cf * scale_m * tmpE(0) * cwtf;
        ky += -cf * scale_m * tmpE(1) * cwtf;
    }

    double dx = RefPartBunch_m->getX0(i);
    double dy = RefPartBunch_m->getY0(i);

    Vector_t KR(0.0, 0.0, 0.0);
    addKR(i, t, KR);
    //FIXME ?? different in bet src
    K += Vector_t(KR(1) * dx + kx, KR(1) * dy + ky, 0.0);
    //
    //K += Vector_t(kx - KR(1) * dx, ky - KR(1) * dy, 0.0);
}


bool RFCavity::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool RFCavity::apply(const Vector_t &R,
                     const Vector_t &P,
                     const double &t,
                     Vector_t &E,
                     Vector_t &B) {
    Vector_t tmpR(R(0), R(1), R(2) - startField_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) >= 0.0 &&
        tmpR(2) < length_m) {
        bool outOfBounds = fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        if (outOfBounds) return true;

        E += (scale_m + scaleError_m) * cos(frequency_m * t + phase_m + phaseError_m) * tmpE;
        B -= (scale_m + scaleError_m) * sin(frequency_m * t + phase_m + phaseError_m) * tmpB;

    }
    return false;
}

bool RFCavity::applyToReferenceParticle(const Vector_t &R,
                                        const Vector_t &P,
                                        const double &t,
                                        Vector_t &E,
                                        Vector_t &B) {

    Vector_t tmpR(R(0), R(1), R(2) - startField_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) >= 0.0 &&
        tmpR(2) < length_m) {
        bool outOfBounds = fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        if (outOfBounds) return true;

        E += scale_m * cos(frequency_m * t + phase_m) * tmpE;
        B -= scale_m * sin(frequency_m * t + phase_m) * tmpB;

    }
    return false;
}

void RFCavity::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    using Physics::two_pi;

    if (bunch == NULL) {
        startField = startField_m;
        endField = endField_m;

        return;
    }

    double rBegin = 0.0, rEnd = 0.0;
    Inform msg("RFCavity ", *gmsg);
    std::stringstream errormsg;
    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);

    fieldmap_m->getFieldDimensions(startField_m, endField_m, rBegin, rEnd);
    if(endField_m > startField_m) {
        msg << level2 << getName() << " using file ";
        fieldmap_m->getInfo(&msg);
        if(std::abs((frequency_m - fieldmap_m->getFrequency()) / frequency_m) > 0.01) {
            errormsg << "FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP '" << filename_m << "';\n"
                     << frequency_m / two_pi * 1e-6 << " MHz <> "
                     << fieldmap_m->getFrequency() / two_pi * 1e-6 << " MHz; TAKE ON THE LATTER";
            std::string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
            ERRORMSG(errormsg_str << "\n" << endl);
            if(Ippl::myNode() == 0) {
                std::ofstream omsg("errormsg.txt", std::ios_base::app);
                omsg << errormsg_str << std::endl;
                omsg.close();
            }
            frequency_m = fieldmap_m->getFrequency();
        }
        length_m = endField_m - startField_m;
        endField = startField + length_m;
    } else {
        endField = startField - 1e-3;
    }
}

// In current version ,this function reads in the cavity voltage profile data from file.
void RFCavity::initialise(PartBunchBase<double, 3> *bunch,
                          std::shared_ptr<AbstractTimeDependence> freq_atd,
                          std::shared_ptr<AbstractTimeDependence> ampl_atd,
                          std::shared_ptr<AbstractTimeDependence> phase_atd) {
    using Physics::pi;

    RefPartBunch_m = bunch;

    /// set the time dependent models
    setAmplitudeModel(ampl_atd);
    setPhaseModel(phase_atd);
    setFrequencyModel(freq_atd);

    ifstream in(filename_m.c_str());
    if(!in.good()) {
        throw GeneralClassicException("RFCavity::initialise",
                                      "failed to open file '" + filename_m + "', please check if it exists");
    }
    *gmsg << "* Read cavity voltage profile data" << endl;

    in >> num_points_m;

    RNormal_m  = std::unique_ptr<double[]>(new double[num_points_m]);
    VrNormal_m = std::unique_ptr<double[]>(new double[num_points_m]);
    DvDr_m     = std::unique_ptr<double[]>(new double[num_points_m]);

    for(int i = 0; i < num_points_m; i++) {
        if(in.eof()) {
            throw GeneralClassicException("RFCavity::initialise",
                                          "not enough data in file '" + filename_m + "', please check the data format");
        }
        in >> RNormal_m[i] >> VrNormal_m[i] >> DvDr_m[i];

        VrNormal_m[i] *= RefPartBunch_m->getQ();
        DvDr_m[i]     *= RefPartBunch_m->getQ();
    }
    sinAngle_m = sin(angle_m / 180.0 * pi);
    cosAngle_m = cos(angle_m / 180.0 * pi);

    if (frequency_name_m != "")
      *gmsg << "* Timedependent frequency model " << frequency_name_m << endl;

    *gmsg << "* Cavity voltage data read successfully!" << endl;
}

void RFCavity::finalise()
{}

bool RFCavity::bends() const {
    return false;
}


void RFCavity::goOnline(const double &) {
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

void RFCavity::setComponentType(std::string name) {
    name = Util::toUpper(name);
    if(name == "STANDING") {
        type_m = SW;
    } else if(name == "SINGLEGAP") {
        type_m = SGSW;
    } else if(name != "") {
        std::stringstream errormsg;
        errormsg << getName() << ": CAVITY TYPE " << name << " DOES NOT EXIST;";
        std::string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
        ERRORMSG(errormsg_str << "\n" << endl);
        if(Ippl::myNode() == 0) {
            ofstream omsg("errormsg.txt", ios_base::app);
            omsg << errormsg_str << endl;
            omsg.close();
        }
        throw GeneralClassicException("RFCavity::setComponentType", errormsg_str);
    } else {
        type_m = SW;
    }

}

string RFCavity::getComponentType()const {
    if(type_m == SGSW)
        return std::string("SINGLEGAP");
    else
        return std::string("STANDING");
}

double RFCavity::getCycFrequency()const {
    return  frequency_m;
}

/**
   \brief used in OPAL-cycl

   Is called from OPAL-cycl and can handle
   time dependent frequency, amplitude and phase

   At the moment (test) only the frequency is time
   dependent

 */
void RFCavity::getMomentaKick(const double normalRadius, double momentum[], const double t, const double dtCorrt, const int PID, const double restMass, const int chargenumber) {
    using Physics::two_pi;
    using Physics::pi;
    using Physics::c;
    double derivate;
    double Voltage;

    double momentum2  = momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2];
    double betgam = sqrt(momentum2);

    double gamma = sqrt(1.0 + momentum2);
    double beta = betgam / gamma;

    Voltage = spline(normalRadius, &derivate) * scale_m * 1.0e6; // V

    double transit_factor = 0.0;
    double Ufactor = 1.0;

    double frequency = frequency_m * frequency_td_m->getValue(t);

    if(gapwidth_m > 0.0) {
        transit_factor = 0.5 * frequency * gapwidth_m * 1.0e-3 / (c * beta);
        Ufactor = sin(transit_factor) / transit_factor;
    }

    Voltage *= Ufactor;

    double nphase = (frequency * (t + dtCorrt) * 1.0e-9) - phi0_m / 180.0 * pi ; // rad/s, ns --> rad
    double dgam = Voltage * cos(nphase) / (restMass);

    double tempdegree = fmod(nphase * 360.0 / two_pi, 360.0);
    if(tempdegree > 270.0) tempdegree -= 360.0;

    gamma += dgam;

    double newmomentum2 = pow(gamma, 2) - 1.0;

    double pr = momentum[0] * cosAngle_m + momentum[1] * sinAngle_m;
    double ptheta = sqrt(newmomentum2 - pow(pr, 2));
    double px = pr * cosAngle_m - ptheta * sinAngle_m ; // x
    double py = pr * sinAngle_m + ptheta * cosAngle_m; // y

    double rotate = -derivate * (scale_m * 1.0e6) / ((rmax_m - rmin_m) / 1000.0) * sin(nphase) / (frequency * two_pi) / (betgam * restMass / c / chargenumber); // radian

    /// B field effects
    momentum[0] =  cos(rotate) * px + sin(rotate) * py;
    momentum[1] = -sin(rotate) * px + cos(rotate) * py;

    if(PID == 0) {

        Inform  m("OPAL", *gmsg, Ippl::myNode());

        m << "* Cavity " << getName() << " Phase= " << tempdegree << " [deg] transit time factor=  " << Ufactor
          << " dE= " << dgam *restMass * 1.0e-6 << " [MeV]"
          << " E_kin= " << (gamma - 1.0)*restMass * 1.0e-6 << " [MeV] Time dep freq = " << frequency_td_m->getValue(t) << endl;
    }

}

/* cubic spline subrutine */
double RFCavity::spline(double z, double *za) {
    double splint;

    // domain-test and handling of case "1-support-point"
    if(num_points_m < 1) {
        throw GeneralClassicException("RFCavity::spline",
                                      "no support points!");
    }
    if(num_points_m == 1) {
        splint = RNormal_m[0];
        *za = 0.0;
        return splint;
    }

    // search the two support-points
    int il, ih;
    il = 0;
    ih = num_points_m - 1;
    while((ih - il) > 1) {
        int i = (int)((il + ih) / 2.0);
        if(z < RNormal_m[i]) {
            ih = i;
        } else if(z > RNormal_m[i]) {
            il = i;
        } else if(z == RNormal_m[i]) {
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

void RFCavity::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


ElementBase::ElementType RFCavity::getType() const {
    return RFCAVITY;
}

double RFCavity::getAutoPhaseEstimateFallback(double E0, double t0, double q, double mass) {
    using Physics::pi;
    const double dt = 1e-13;
    const double p0 = Util::getP(E0, mass);
    const double origPhase =getPhasem();

    double dphi = pi / 18;

    double phi = 0.0;
    setPhasem(phi);
    std::pair<double, double> ret = trackOnAxisParticle(E0 / mass, t0, dt, q, mass);
    double phimax = 0.0;
    double Emax = Util::getEnergy(Vector_t(0.0, 0.0, ret.first), mass);
    phi += dphi;

    for (unsigned int j = 0; j < 2; ++ j) {
        for (unsigned int i = 0; i < 36; ++ i, phi += dphi) {
            setPhasem(phi);
            ret = trackOnAxisParticle(p0, t0, dt, q, mass);
            double Ekin = Util::getEnergy(Vector_t(0.0, 0.0, ret.first), mass);
            if (Ekin > Emax) {
                Emax = Ekin;
                phimax = phi;
            }
        }

        phi = phimax - dphi;
        dphi = dphi / 17.5;
    }

    phimax = phimax - floor(phimax / Physics::two_pi + 0.5) * Physics::two_pi;
    phimax = fmod(phimax, Physics::two_pi);

    const int prevPrecision = Ippl::Info->precision(8);
    INFOMSG(level2
            << "estimated phase= " << phimax << " rad = "
            << phimax * Physics::rad2deg << " deg \n"
            << "Ekin= " << Emax << " MeV" << setprecision(prevPrecision) << "\n" << endl);

    setPhasem(origPhase);
    return phimax;
}

double RFCavity::getAutoPhaseEstimate(const double &E0, const double &t0, const double &q, const double &mass) {
    vector<double> t, E, t2, E2;
    std::vector< double > F;
    std::vector< std::pair< double, double > > G;
    gsl_spline *onAxisInterpolants;
    gsl_interp_accel *onAxisAccel;

    unsigned int N;
    double A, B;
    double phi = 0.0, tmp_phi, dphi = 0.5 * Physics::pi / 180.;
    double dz = 1.0, length = 0.0;
    fieldmap_m->getOnaxisEz(G);
    double begin = (G.front()).first;
    double end = (G.back()).first;
    std::unique_ptr<double[]> zvals(new double[G.size()]);
    std::unique_ptr<double[]> onAxisField(new double[G.size()]);

    for(size_t j = 0; j < G.size(); ++ j) {
        zvals[j] = G[j].first;
        onAxisField[j] = G[j].second;
    }
    onAxisInterpolants = gsl_spline_alloc(gsl_interp_cspline, G.size());
    onAxisAccel = gsl_interp_accel_alloc();
    gsl_spline_init(onAxisInterpolants, zvals.get(), onAxisField.get(), G.size());

    length = end - begin;
    dz = length / G.size();

    G.clear();

    N = (int)floor(length / dz + 1);
    dz = length / N;

    F.resize(N);
    double z = begin;
    for(size_t j = 0; j < N; ++ j, z += dz) {
        F[j] = gsl_spline_eval(onAxisInterpolants, z, onAxisAccel);
    }
    gsl_spline_free(onAxisInterpolants);
    gsl_interp_accel_free(onAxisAccel);

    t.resize(N, t0);
    t2.resize(N, t0);
    E.resize(N, E0);
    E2.resize(N, E0);

    z = begin + dz;
    for(unsigned int i = 1; i < N; ++ i, z += dz) {
        E[i] = E[i - 1] + dz * scale_m / mass;
        E2[i] = E[i];
    }

    for(int iter = 0; iter < 10; ++ iter) {
        A = B = 0.0;
        for(unsigned int i = 1; i < N; ++ i) {
            t[i] = t[i - 1] + getdT(i, E, dz, mass);
            t2[i] = t2[i - 1] + getdT(i, E2, dz, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, t, dz, frequency_m, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, t, dz, frequency_m, F);
        }

        if(std::abs(B) > 0.0000001) {
            tmp_phi = atan(A / B);
        } else {
            tmp_phi = Physics::pi / 2;
        }
        if(q * (A * sin(tmp_phi) + B * cos(tmp_phi)) < 0) {
            tmp_phi += Physics::pi;
        }

        if(std::abs(phi - tmp_phi) < frequency_m * (t[N - 1] - t[0]) / (10 * N)) {
            for(unsigned int i = 1; i < N; ++ i) {
                E[i] = E[i - 1];
                E[i] += q * scale_m * getdE(i, t, dz, phi, frequency_m, F) ;
            }
            const int prevPrecision = Ippl::Info->precision(8);
            INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad = "
                    << tmp_phi * Physics::rad2deg << " deg \n"
                    << "Ekin= " << E[N - 1] << " MeV" << setprecision(prevPrecision) << "\n" << endl);

            return tmp_phi;
        }
        phi = tmp_phi - floor(tmp_phi / Physics::two_pi + 0.5) * Physics::two_pi;

        for(unsigned int i = 1; i < N; ++ i) {
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

            E[i] = E[i - 1];
            E2[i] = E2[i - 1];
            E[i] += q * scale_m * getdE(i, t, dz, phi, frequency_m, F) ;
            E2[i] += q * scale_m * getdE(i, t2, dz, phi + dphi, frequency_m, F);
        }

        double cosine_part = 0.0, sine_part = 0.0;
        double p0 = sqrt((E0 / mass + 1) * (E0 / mass + 1) - 1);
        cosine_part += scale_m * cos(frequency_m * t0) * F[0];
        sine_part += scale_m * sin(frequency_m * t0) * F[0];

        double totalEz0 = cos(phi) * cosine_part - sin(phi) * sine_part;

        if(p0 + q * totalEz0 * (t[1] - t[0]) * Physics::c / mass < 0) {
            // make totalEz0 = 0
            tmp_phi = atan(cosine_part / sine_part);
            if(abs(tmp_phi - phi) > Physics::pi) {
                phi = tmp_phi + Physics::pi;
            } else {
                phi = tmp_phi;
            }
        }
    }

    const int prevPrecision = Ippl::Info->precision(8);
    INFOMSG(level2
            << "estimated phase= " << tmp_phi << " rad = "
            << tmp_phi * Physics::rad2deg << " deg \n"
            << "Ekin= " << E[N - 1] << " MeV" << setprecision(prevPrecision) << "\n" << endl);

    return phi;
}

pair<double, double> RFCavity::trackOnAxisParticle(const double &p0,
                                                   const double &t0,
                                                   const double &dt,
                                                   const double &q,
                                                   const double &mass,
                                                   std::ofstream *out) {
    Vector_t p(0, 0, p0);
    double t = t0;
    BorisPusher integrator(*RefPartBunch_m->getReference());
    const double cdt = Physics::c * dt;
    const double zbegin = startField_m;
    const double zend = length_m + startField_m;

    Vector_t z(0.0, 0.0, zbegin);
    double dz = 0.5 * p(2) / sqrt(1.0 + dot(p, p)) * cdt;
    Vector_t Ef(0.0), Bf(0.0);

    if (out) *out << std::setw(18) << z[2]
                  << std::setw(18) << Util::getEnergy(p, mass)
                  << std::endl;
    while(z(2) + dz < zend && z(2) + dz > zbegin) {
        z /= cdt;
        integrator.push(z, p, dt);
        z *= cdt;

        if(z(2) >= zbegin && z(2) <= zend) {
            Ef = 0.0;
            Bf = 0.0;
            applyToReferenceParticle(z, p, t + 0.5 * dt, Ef, Bf);
        }
        integrator.kick(z, p, Ef, Bf, dt);

        dz = 0.5 * p(2) / sqrt(1.0 + dot(p, p)) * cdt;
        z /= cdt;
        integrator.push(z, p, dt);
        z *= cdt;
        t += dt;

        if (out) *out << std::setw(18) << z[2]
                      << std::setw(18) << Util::getEnergy(p, mass)
                      << std::endl;
    }

    const double beta = sqrt(1. - 1 / (dot(p, p) + 1.));
    const double tErr  = (z(2) - zend) / (Physics::c * beta);

    return pair<double, double>(p(2), t - tErr);
}

bool RFCavity::isInside(const Vector_t &r) const {
    if (isInsideTransverse(r)) {
        return fieldmap_m->isInside(r);
    }

    return false;
}

double RFCavity::getElementLength() const {
    double length = ElementBase::getElementLength();
    if (length < 1e-10 && fieldmap_m != NULL) {
        double start, end, tmp;
        fieldmap_m->getFieldDimensions(start, end, tmp, tmp);
        length = end - start;
    }

    return length;
}

void RFCavity::getElementDimensions(double &begin,
                                    double &end) const {
    double tmp;
    fieldmap_m->getFieldDimensions(begin, end, tmp, tmp);
}