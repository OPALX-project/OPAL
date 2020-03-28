// ------------------------------------------------------------------------
// $RCSfile: TravelingWave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TravelingWave
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

#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

#include <cmath>
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class TravelingWave
// ------------------------------------------------------------------------

TravelingWave::TravelingWave():
    TravelingWave("")
{}


TravelingWave::TravelingWave(const TravelingWave &right):
    RFCavity(right),
    CoreFieldmap_m(NULL),
    scaleCore_m(right.scaleCore_m),
    scaleCoreError_m(right.scaleCoreError_m),
    phaseCore1_m(right.phaseCore1_m),
    phaseCore2_m(right.phaseCore2_m),
    phaseExit_m(right.phaseExit_m),
    length_m(right.length_m),
    startCoreField_m(right.startCoreField_m),
    startExitField_m(right.startExitField_m),
    mappedStartExitField_m(right.mappedStartExitField_m),
    PeriodLength_m(right.PeriodLength_m),
    NumCells_m(right.NumCells_m),
    CellLength_m(right.CellLength_m),
    Mode_m(right.Mode_m),
    fast_m(right.fast_m),
    autophaseVeto_m(right.autophaseVeto_m),
    designEnergy_m(right.designEnergy_m)
{}


TravelingWave::TravelingWave(const std::string &name):
    RFCavity(name),
    CoreFieldmap_m(NULL),
    scaleCore_m(1.0),
    scaleCoreError_m(0.0),
    phaseCore1_m(0.0),
    phaseCore2_m(0.0),
    phaseExit_m(0.0),
    length_m(0.0),
    startCoreField_m(0.0),
    startExitField_m(0.0),
    mappedStartExitField_m(0.0),
    PeriodLength_m(0.0),
    NumCells_m(1),
    CellLength_m(0.0),
    Mode_m(1),
    fast_m(true),
    autophaseVeto_m(false),
    designEnergy_m(-1.0)
{}


TravelingWave::~TravelingWave() {
    // Fieldmap::deleteFieldmap(filename_m);
    //   Fieldmap::deleteFieldmap(EntryFilename_m);
    //   Fieldmap::deleteFieldmap(ExitFilename_m);
}


void TravelingWave::accept(BeamlineVisitor &visitor) const {
    visitor.visitTravelingWave(*this);
}

/**
 * ENVELOPE COMPONENT for radial focussing of the beam
 * Calculates the transverse envelope component for the traveling wave
 * element and adds it to the K vector
*/
void TravelingWave::addKR(int i, double t, Vector_t &K) {
    Inform msg("TravelingWave::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0), tmpB_diff(0.0, 0.0, 0.0);
    Vector_t tmpA0(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i) - 0.5 * PeriodLength_m);

    double g = RefPartBunch_m->getGamma(i);
    double b = RefPartBunch_m->getBeta(i);
    DiffDirection zdir(DZ);
    double wtf = 0.0;
    double k = 0.0;

    if(tmpA0(2) < startCoreField_m) {

        wtf = frequency_m * t + phase_m;
        CoreFieldmap_m->getFieldstrength(tmpA0, tmpE, tmpB);
        CoreFieldmap_m->getFieldDerivative(tmpA0, tmpE_diff, tmpB_diff, zdir);
        k = scale_m * (tmpE_diff(2) * std::cos(wtf) - b * frequency_m * tmpE(2) * std::sin(wtf) / Physics::c);

    } else if(tmpA0(2) < startExitField_m) {

        wtf = frequency_m * t + phaseCore1_m;
        tmpA0(2) -= startCoreField_m;
        const double z = tmpA0(2);
        tmpA0(2) = tmpA0(2) - PeriodLength_m * std::floor(tmpA0(2) / PeriodLength_m);
        tmpA0(2) += startCoreField_m;
        tmpE = 0.0;
        tmpB = 0.0;
        tmpE_diff = 0.0;
        tmpB_diff = 0.0;

        CoreFieldmap_m->getFieldstrength(tmpA0, tmpE, tmpB);
        CoreFieldmap_m->getFieldDerivative(tmpA0, tmpE_diff, tmpB_diff, zdir);
        k = scaleCore_m * (tmpE_diff(2) * std::cos(wtf) - b * frequency_m * tmpE(2) * std::sin(wtf) / Physics::c);

        wtf = frequency_m * t + phaseCore2_m;
        tmpA0(2) = z + CellLength_m;
        tmpA0(2) = tmpA0(2) - PeriodLength_m * std::floor(tmpA0(2) / PeriodLength_m);
        tmpA0(2) += startCoreField_m;
        tmpE = 0.0;
        tmpB = 0.0;
        tmpE_diff = 0.0;
        tmpB_diff = 0.0;

        CoreFieldmap_m->getFieldstrength(tmpA0, tmpE, tmpB);
        CoreFieldmap_m->getFieldDerivative(tmpA0, tmpE_diff, tmpB_diff, zdir);
        k += scaleCore_m * (tmpE_diff(2) * std::cos(wtf) - b * frequency_m * tmpE(2) * std::sin(wtf) / Physics::c);

    } else {

        wtf  = frequency_m * t + phaseExit_m;
        tmpA0(2) -= mappedStartExitField_m;
        tmpE = 0.0;
        tmpB = 0.0;
        tmpE_diff = 0.0;
        tmpB_diff = 0.0;

        CoreFieldmap_m->getFieldstrength(tmpA0, tmpE, tmpB);
        CoreFieldmap_m->getFieldDerivative(tmpA0, tmpE_diff, tmpB_diff, zdir);
        k = scale_m * (tmpE_diff(2) * std::cos(wtf) - b * frequency_m * tmpE(2) * std::sin(wtf) / Physics::c);

    }

    //XXX: we have to compensate for bet behaviour by using -q_e
    k *= (-Physics::q_e / (2.0 * g * Physics::EMASS));
    K += Vector_t(k, k, 0.0);
}


void TravelingWave::addKT(int i, double t, Vector_t &K) {
    // get K from radial function, solves the lastk problem from BET
    Vector_t tempK(0.0, 0.0, 0.0);
    addKR(i, t, tempK);

    // x and y component are identical and equal to k
    double oldK = tempK(0);

    //get coordinates
    double dx = RefPartBunch_m->getX0(i);
    double dy = RefPartBunch_m->getY0(i);

    K += Vector_t(oldK * dx, oldK * dy, 0.0);
}

/*
    Function to obtain z-k-component, copied from BET:
   Field Kloc;

   if (isActive(z)){
       double k,kx,ky;

       if (kLast == MAXFLOAT) {
           Kloc = getK(g,b,z,t,m);  // dummy call to set kLast
       }
       k   = kLast;
       kx  = 0.0; ky = 0.0;
       if (Cxy < 1.0) {
           Field Eloc = getE(z,t);
           double
               cf = -e0/(g*m);

           kx = cf*Eloc.x;
           ky = cf*Eloc.y;
       }
       Kloc = Field(k*Cxy*(x0-x) + kx,
               k*Cxy*(y0-y) + ky,
               0.0);
       if (Ndamp > 0) {
           Kloc *= fDamp(z - z0);
       }
   }

   return Kloc;
*/


bool TravelingWave::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool TravelingWave::apply(const Vector_t &R, const Vector_t &/*P*/, const double &t, Vector_t &E, Vector_t &B) {
    if (R(2) < -0.5 * PeriodLength_m || R(2) + 0.5 * PeriodLength_m >= length_m) return false;

    Vector_t tmpR = Vector_t(R(0), R(1), R(2) + 0.5 * PeriodLength_m);
    double tmpcos, tmpsin;
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) < startCoreField_m) {
        if (!CoreFieldmap_m->isInside(tmpR)) return true;

        tmpcos =  (scale_m + scaleError_m) * std::cos(frequency_m * t + phase_m + phaseError_m);
        tmpsin = -(scale_m + scaleError_m) * std::sin(frequency_m * t + phase_m + phaseError_m);

    } else if (tmpR(2) < startExitField_m) {
        Vector_t tmpE2(0.0, 0.0, 0.0), tmpB2(0.0, 0.0, 0.0);
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - PeriodLength_m * std::floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;
        if (!CoreFieldmap_m->isInside(tmpR)) return true;

        tmpcos =  (scaleCore_m + scaleCoreError_m) * std::cos(frequency_m * t + phaseCore1_m + phaseError_m);
        tmpsin = -(scaleCore_m + scaleCoreError_m) * std::sin(frequency_m * t + phaseCore1_m + phaseError_m);
        CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + CellLength_m;
        tmpR(2) = tmpR(2) - PeriodLength_m * std::floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  (scaleCore_m + scaleCoreError_m) * std::cos(frequency_m * t + phaseCore2_m + phaseError_m);
        tmpsin = -(scaleCore_m + scaleCoreError_m) * std::sin(frequency_m * t + phaseCore2_m + phaseError_m);

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!CoreFieldmap_m->isInside(tmpR)) return true;
        tmpcos =  (scale_m + scaleError_m) * std::cos(frequency_m * t + phaseExit_m + phaseError_m);
        tmpsin = -(scale_m + scaleError_m) * std::sin(frequency_m * t + phaseExit_m + phaseError_m);
    }

    CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return false;
}

bool TravelingWave::applyToReferenceParticle(const Vector_t &R, const Vector_t &/*P*/, const double &t, Vector_t &E, Vector_t &B) {

    if (R(2) < -0.5 * PeriodLength_m || R(2) + 0.5 * PeriodLength_m >= length_m) return false;

    Vector_t tmpR = Vector_t(R(0), R(1), R(2) + 0.5 * PeriodLength_m);
    double tmpcos, tmpsin;
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) < startCoreField_m) {
        if (!CoreFieldmap_m->isInside(tmpR)) return true;
        tmpcos =  scale_m * std::cos(frequency_m * t + phase_m);
        tmpsin = -scale_m * std::sin(frequency_m * t + phase_m);

    } else if (tmpR(2) < startExitField_m) {
        Vector_t tmpE2(0.0, 0.0, 0.0), tmpB2(0.0, 0.0, 0.0);
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - PeriodLength_m * std::floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;
        if (!CoreFieldmap_m->isInside(tmpR)) return true;

        tmpcos =  scaleCore_m * std::cos(frequency_m * t + phaseCore1_m);
        tmpsin = -scaleCore_m * std::sin(frequency_m * t + phaseCore1_m);
        CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + CellLength_m;
        tmpR(2) = tmpR(2) - PeriodLength_m * std::floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  scaleCore_m * std::cos(frequency_m * t + phaseCore2_m);
        tmpsin = -scaleCore_m * std::sin(frequency_m * t + phaseCore2_m);

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!CoreFieldmap_m->isInside(tmpR)) return true;

        tmpcos =  scale_m * std::cos(frequency_m * t + phaseExit_m);
        tmpsin = -scale_m * std::sin(frequency_m * t + phaseExit_m);
    }

    CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return false;
}

void TravelingWave::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    using Physics::pi;
    using Physics::two_pi;

    if (bunch == NULL) {
        startField = - 0.5 * PeriodLength_m;
        endField = startExitField_m;
        return;
    }

    Inform msg("TravelingWave ", *gmsg);
    std::stringstream errormsg;

    RefPartBunch_m = bunch;

    CoreFieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    if(CoreFieldmap_m != NULL) {
        double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
        CoreFieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

        if(zEnd > zBegin) {
            msg << level2 << getName() << " using file ";
            CoreFieldmap_m->getInfo(&msg);
            if(std::abs((frequency_m - CoreFieldmap_m->getFrequency()) / frequency_m) > 0.01) {
                errormsg << "FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP '" <<  filename_m + "';\n"
                         << frequency_m / two_pi * 1e-6 << " MHz <> "
                         << CoreFieldmap_m->getFrequency() / two_pi * 1e-6 << " MHz; TAKE ON THE LATTER\n";
                std::string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
                *gmsg << errormsg_str << endl;
                ERRORMSG(errormsg_str << "\n" << endl);
                if(Ippl::myNode() == 0) {
                    std::ofstream omsg("errormsg.txt", std::ios_base::app);
                    omsg << errormsg_str << std::endl;
                    omsg.close();
                }
                frequency_m = CoreFieldmap_m->getFrequency();
            }

            PeriodLength_m = (zEnd - zBegin) / 2.0;
            CellLength_m = PeriodLength_m * Mode_m;

            startCoreField_m = PeriodLength_m / 2.0;
            startExitField_m = startCoreField_m + (NumCells_m - 1) * CellLength_m;
            mappedStartExitField_m = startExitField_m - 3.0 * PeriodLength_m / 2.0;

            startField = -PeriodLength_m / 2.0;
            endField = startField + startExitField_m + PeriodLength_m / 2.0;
            length_m = endField - startField;

            scaleCore_m = scale_m / std::sin(2.0 * pi * Mode_m);
            scaleCoreError_m = scaleError_m / std::sin(2.0 * pi * Mode_m);
            phaseCore1_m = phase_m + pi * Mode_m / 2.0;
            phaseCore2_m = phase_m + pi * Mode_m * 1.5;
            phaseExit_m = phase_m - 2.0 * pi * ((NumCells_m - 1) * Mode_m - std::floor((NumCells_m - 1) * Mode_m));

        } else {
            endField = startField - 1e-3;
        }
    } else {
        endField = startField - 1e-3;
    }
}

void TravelingWave::finalise()
{}

bool TravelingWave::bends() const {
    return false;
}


void TravelingWave::goOnline(const double &) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void TravelingWave::goOffline() {
    Fieldmap::freeMap(filename_m);
}

void TravelingWave::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = -0.5 * PeriodLength_m;
    zEnd = zBegin + length_m;
}


double TravelingWave::getElementLength() const {
    return length_m;
}

void TravelingWave::getElementDimensions(double &begin,
                                         double &end) const {
    begin = -0.5 * PeriodLength_m;
    end = begin + length_m;
}

ElementBase::ElementType TravelingWave::getType() const {
    return TRAVELINGWAVE;
}

double TravelingWave::getAutoPhaseEstimate(const double &E0, const double &t0, const double &q, const double &mass) {
    std::vector<double> t, E, t2, E2;
    std::vector<std::pair<double, double> > F;
    double Dz;
    int N1, N2, N3, N4;
    double A, B;
    double phi = 0.0, tmp_phi, dphi = 0.5 * Physics::pi / 180.;
    double phaseC1 = phaseCore1_m - phase_m;
    double phaseC2 = phaseCore2_m - phase_m;
    double phaseE = phaseExit_m - phase_m;

    CoreFieldmap_m->getOnaxisEz(F);
    if(F.size() == 0) return 0.0;

    N1 = static_cast<int>(std::floor(F.size() / 4.)) + 1;
    N2 = F.size() - 2 * N1 + 1;
    N3 = 2 * N1 + static_cast<int>(std::floor((NumCells_m - 1) * N2 * Mode_m)) - 1;
    N4 = static_cast<int>(std::round(N2 * Mode_m));
    Dz = F[N1 + N2].first - F[N1].first;

    t.resize(N3, t0);
    t2.resize(N3, t0);
    E.resize(N3, E0);
    E2.resize(N3, E0);
    for(int i = 1; i < N1; ++ i) {
        E[i] = E0 + (F[i].first + F[i - 1].first) / 2. * scale_m / mass;
        E2[i] = E[i];
    }
    for(int i = N1; i < N3 - N1 + 1; ++ i) {
        int I = (i - N1) % N2 + N1;
        double Z = (F[I].first + F[I - 1].first) / 2. + std::floor((i - N1) / N2) * Dz;
        E[i] = E0 + Z * scaleCore_m / mass;
        E2[i] = E[i];
    }
    for(int i = N3 - N1 + 1; i < N3; ++ i) {
        int I = i - N3 - 1 + 2 * N1 + N2;
        double Z = (F[I].first + F[I - 1].first) / 2. + ((NumCells_m - 1) * Mode_m - 1) * Dz;
        E[i] = E0 + Z * scale_m / mass;
        E2[i] = E[i];
    }

    for(int iter = 0; iter < 10; ++ iter) {
        A = B = 0.0;
        for(int i = 1; i < N1; ++ i) {
            t[i] = t[i - 1] + getdT(i, i, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, i, E2, F, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, i, t, 0.0, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, i, t, 0.0, F);
        }
        for(int i = N1; i < N3 - N1 + 1; ++ i) {
            int I = (i - N1) % N2 + N1;
            int J = (i - N1 + N4) % N2 + N1;
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            A += scaleCore_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * (getdA(i, I, t, phaseC1, F) + getdA(i, J, t, phaseC2, F));
            B += scaleCore_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * (getdB(i, I, t, phaseC1, F) + getdB(i, J, t, phaseC2, F));
        }
        for(int i = N3 - N1 + 1; i < N3; ++ i) {
            int I = i - N3 - 1 + 2 * N1 + N2;
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, I, t, phaseE, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, I, t, phaseE, F);
        }

        if(std::abs(B) > 0.0000001) {
            tmp_phi = std::atan(A / B);
        } else {
            tmp_phi = Physics::pi / 2;
        }
        if(q * (A * std::sin(tmp_phi) + B * std::cos(tmp_phi)) < 0) {
            tmp_phi += Physics::pi;
        }

        if(std::abs(phi - tmp_phi) < frequency_m * (t[N3 - 1] - t[0]) / N3) {
            for(int i = 1; i < N1; ++ i) {
                E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            }
            for(int i = N1; i < N3 - N1 + 1; ++ i) {
                int I = (i - N1) % N2 + N1;
                int J = (i - N1 + N4) % N2 + N1;
                E[i] = E[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            }
            for(int i = N3 - N1 + 1; i < N3; ++ i) {
                int I = i - N3 - 1 + 2 * N1 + N2;
                E[i] = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            }

            const int prevPrecision = Ippl::Info->precision(8);
            INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad = "
                    << tmp_phi * Physics::rad2deg << " deg,\n"
                    << "Ekin= " << E[N3 - 1] << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);
            return tmp_phi;
        }
        phi = tmp_phi - std::round(tmp_phi / Physics::two_pi) * Physics::two_pi;


        for(int i = 1; i < N1; ++ i) {
            E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, i, t, phi + dphi, F); // should I use here t or t2?
            t[i] = t[i - 1] + getdT(i, i, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, i, E2, F, mass);
            E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, i, t2, phi + dphi, F);
        }
        for(int i = N1; i < N3 - N1 + 1; ++ i) {
            int I = (i - N1) % N2 + N1;
            int J = (i - N1 + N4) % N2 + N1;
            E[i] = E[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            E2[i] = E2[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1 + dphi, F) + getdE(i, J, t, phi + phaseC2 + dphi, F)); //concerning t: see above
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            E[i] = E[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            E2[i] = E2[i - 1] + q * scaleCore_m * (getdE(i, I, t2, phi + phaseC1 + dphi, F) + getdE(i, J, t2, phi + phaseC2 + dphi, F));
        }
        for(int i = N3 - N1 + 1; i < N3; ++ i) {
            int I = i - N3 - 1 + 2 * N1 + N2;
            E[i] = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE + dphi, F); //concerning t: see above
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            E[i] = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, I, t2, phi + phaseE + dphi, F);
        }
        //             msg << ", Ekin= " << E[N3-1] << " MeV" << endl;
    }


    const int prevPrecision = Ippl::Info->precision(8);
    INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad = "
            << tmp_phi * Physics::rad2deg << " deg,\n"
            << "Ekin= " << E[N3 - 1] << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);

    return phi;
}

std::pair<double, double> TravelingWave::trackOnAxisParticle(const double &p0,
                                                             const double &t0,
                                                             const double &dt,
                                                             const double &q,
                                                             const double &mass,
                                                             std::ofstream *out) {
    double phase = frequency_m * t0 + phase_m;
    double p = p0;
    double t = t0;
    double cdt = Physics::c * dt;
    double dphi = frequency_m * dt;
    std::vector<std::pair<double, double> > F;
    CoreFieldmap_m->getOnaxisEz(F);

    double *zvals = new double[F.size()];
    double *onAxisField = new double[F.size()];
    double Ezmax = 0.0;
    for(unsigned int i = 0; i < F.size(); ++i) {
        zvals[i] = F[i].first;
        onAxisField[i] = F[i].second;
        if(std::abs(onAxisField[i]) > Ezmax) Ezmax = std::abs(onAxisField[i]);
    }
    //    Ezmax /= 1e6;
    double z = zvals[0];
    double zbegin = z;

    gsl_spline *onAxisInterpolants = gsl_spline_alloc(gsl_interp_cspline, F.size());
    gsl_interp_accel *onAxisAccel = gsl_interp_accel_alloc();
    gsl_spline_init(onAxisInterpolants, zvals, onAxisField, F.size());

    delete[] zvals;
    delete[] onAxisField;

    if (out) *out << std::setw(18) << z
                  << std::setw(18) << Util::getEnergy(p, mass)
                  << std::endl;
    double dz = 0.5 * p / std::sqrt(1.0 + p * p) * cdt;
    while(z + dz < startCoreField_m + zbegin) {
        z += dz;
        double ez = scale_m / Ezmax * std::cos(phase) * gsl_spline_eval(onAxisInterpolants, z, onAxisAccel);
        p += ez * q * cdt / mass;
        dz = 0.5 * p / std::sqrt(1.0 + p * p) * cdt;
        z += 0.5 * p / std::sqrt(1.0 + p * p) * cdt;
        phase += dphi;
        t += dt;
    }
    double phase2 = phase - phase_m + phaseCore2_m;
    phase = phase - phase_m + phaseCore1_m;
    while(z + dz < startExitField_m + zbegin) {
        z += dz;

        double tmpz = z - (startCoreField_m + zbegin);
        tmpz -= PeriodLength_m * std::floor(tmpz / PeriodLength_m);
        tmpz += startCoreField_m + zbegin;

        double ez = scaleCore_m / Ezmax * std::cos(phase) * gsl_spline_eval(onAxisInterpolants, tmpz, onAxisAccel);

        tmpz = z - (startCoreField_m + zbegin) + CellLength_m;
        tmpz -= PeriodLength_m * std::floor(tmpz / PeriodLength_m);
        tmpz += startCoreField_m + zbegin;

        ez += scaleCore_m / Ezmax * std::cos(phase2) * gsl_spline_eval(onAxisInterpolants, tmpz, onAxisAccel);

        p += ez * q * cdt / mass;
        dz = 0.5 * p / std::sqrt(1.0 + p * p) * cdt;
        z += dz;
        phase += dphi;
        phase2 += dphi;
        t += dt;

        if (out) *out << std::setw(18) << z
                      << std::setw(18) << Util::getEnergy(p, mass)
                      << std::endl;
    }
    phase = phase - phaseCore1_m + phaseExit_m;
    while(z + dz < startExitField_m + 0.5 * PeriodLength_m + zbegin) {
        z += dz;
        double tmpz = z - (startExitField_m + zbegin);
        tmpz -= PeriodLength_m * std::floor(tmpz / PeriodLength_m);
        tmpz += 3.0 * PeriodLength_m / 2.0;
        double ez = scale_m / Ezmax * std::cos(phase) * gsl_spline_eval(onAxisInterpolants, tmpz, onAxisAccel);
        p += ez * q * cdt / mass;
        dz = 0.5 * p / std::sqrt(1.0 + p * p) * cdt;
        z += dz;
        phase += dphi;
        t += dt;
    }

    gsl_spline_free(onAxisInterpolants);
    gsl_interp_accel_free(onAxisAccel);

    const double beta = std::sqrt(1. - 1 / (p * p + 1.));
    const double tErr  = (z - (startExitField_m + 0.5 * PeriodLength_m + zbegin)) / (Physics::c * beta);

    return std::pair<double, double>(p, t - tErr);
}

bool TravelingWave::isInside(const Vector_t &r) const {
    return (isInsideTransverse(r) && r(2) >= -0.5 * PeriodLength_m && r(2) < startExitField_m);
}
