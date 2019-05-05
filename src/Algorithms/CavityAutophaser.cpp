#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"

#include <fstream>
#include <iostream>

extern Inform *gmsg;

CavityAutophaser::CavityAutophaser(const PartData &ref,
                                   std::shared_ptr<Component> cavity):
    itsReference_m(ref),
    itsCavity_m(cavity)
{
    double zbegin = 0.0, zend = 0.0;
    PartBunchBase<double, 3> *fakeBunch = NULL;
    cavity->initialise(fakeBunch, zbegin, zend);
    initialR_m = Vector_t(0, 0, zbegin);
}

CavityAutophaser::~CavityAutophaser() {

}

double CavityAutophaser::getPhaseAtMaxEnergy(const Vector_t &R,
                                             const Vector_t &P,
                                             double t,
                                             double dt) {
    if(!(itsCavity_m->getType() == ElementBase::TRAVELINGWAVE ||
         itsCavity_m->getType() == ElementBase::RFCAVITY)) {
        throw OpalException("CavityAutophaser::getPhaseAtMaxEnergy()",
                            "given element is not a cavity");
    }

    initialP_m = Vector_t(0, 0, euclidean_norm(P));

    RFCavity *element     = static_cast<RFCavity *>(itsCavity_m.get());
    bool apVeto           = element->getAutophaseVeto();
    bool isDCGun          = false;
    double originalPhase  = element->getPhasem();
    double tErr           = (initialR_m(2) - R(2)) * sqrt(dot(P,P) + 1.0) / (P(2) * Physics::c);
    double optimizedPhase = 0.0;
    double finalEnergy    = 0.0;
    double newPhase       = 0.0;
    double amplitude      = element->getAmplitudem();
    double basePhase      = std::fmod(element->getFrequencym() * (t + tErr), Physics::two_pi);
    double frequency      = element->getFrequencym();

    if ((!apVeto) && frequency <= (1.0 + 1e-6) * Physics::two_pi) { // DC gun
        optimizedPhase = (amplitude * itsReference_m.getQ() > 0.0? 0.0: Physics::pi);
        element->setPhasem(optimizedPhase + originalPhase);
        element->setAutophaseVeto();

        originalPhase += optimizedPhase;
        OpalData::getInstance()->setMaxPhase(itsCavity_m->getName(), originalPhase);

        apVeto = true;
        isDCGun = true;
    }

    std::stringstream ss;
    for (char c: itsCavity_m->getName()) {
        ss << std::setw(2) << std::left << c;
    }
    INFOMSG(level1 << "\n* ************* "
                   << std::left << std::setw(68) << std::setfill('*') << ss.str()
                   << std::setfill(' ') << endl);
    if (!apVeto) {
        double initialEnergy = Util::getEnergy(P, itsReference_m.getM()) * 1e-6;
        double AstraPhase    = 0.0;
        double designEnergy  = element->getDesignEnergy();

        if (amplitude < 0.0) {
            amplitude = -amplitude;
            element->setAmplitudem(amplitude);
        }

        double initialPhase  = guessCavityPhase(t + tErr);

        if (amplitude == 0.0 && designEnergy <= 0.0) {
            throw OpalException("CavityAutophaser::getPhaseAtMaxEnergy()",
                                "neither amplitude or design energy given to cavity " + element->getName());
        }

        if (designEnergy > 0.0) {
            const double length = itsCavity_m->getElementLength();
            if (length <= 0.0) {
                throw OpalException("CavityAutophaser::getPhaseAtMaxEnergy()",
                                    "length of cavity " + element->getName() + " is zero");
            }

            amplitude = 2 * (designEnergy - initialEnergy) / (std::abs(itsReference_m.getQ()) * length);

            element->setAmplitudem(amplitude);

            int count = 0;
            while (count < 1000) {
                initialPhase = guessCavityPhase(t + tErr);
                auto status = optimizeCavityPhase(initialPhase, t + tErr, dt);

                optimizedPhase = status.first;
                finalEnergy = status.second;

                if (std::abs(designEnergy - finalEnergy) < 1e-7) break;

                amplitude *= std::abs(designEnergy / finalEnergy);
                element->setAmplitudem(amplitude);
                initialPhase = optimizedPhase;

                ++ count;
            }
        }
        auto status = optimizeCavityPhase(initialPhase, t + tErr, dt);

        optimizedPhase = status.first;
        finalEnergy = status.second;

        AstraPhase = std::fmod(optimizedPhase + Physics::pi / 2 + Physics::two_pi, Physics::two_pi);
        newPhase = std::fmod(originalPhase + optimizedPhase + Physics::two_pi, Physics::two_pi);
        element->setPhasem(newPhase);
        element->setAutophaseVeto();
        OpalData::getInstance()->setMaxPhase(itsCavity_m->getName(), newPhase);

        newPhase = std::fmod(newPhase + basePhase, Physics::two_pi);

        auto opal = OpalData::getInstance();
        if (!opal->isOptimizerRun()) {
            std::ofstream out("data/" + itsCavity_m->getName() + "_AP.dat");
            track(initialR_m, initialP_m, t + tErr, dt, newPhase, &out);
            out.close();
        } else {
            track(initialR_m, initialP_m, t + tErr, dt, newPhase, NULL);
        }

        INFOMSG(level1 << std::fixed << std::setprecision(4)
                << itsCavity_m->getName() << "_phi = "  << newPhase * Physics::rad2deg <<  " [deg], "
                << "corresp. in Astra = " << AstraPhase * Physics::rad2deg << " [deg],\n"
                << "E = " << finalEnergy << " [MeV], " << "phi_nom = " << originalPhase * Physics::rad2deg << " [deg]\n"
                << "Ez_0 = " << amplitude << " [MV/m]" << "\n"
                << "time = " << (t + tErr) * 1e9 << " [ns], dt = " << dt * 1e12 << " [ps]" << endl);

    } else {
        auto status = optimizeCavityPhase(originalPhase, t + tErr, dt);

        finalEnergy = status.second;

        originalPhase = std::fmod(originalPhase, Physics::two_pi);
        double AstraPhase = std::fmod(optimizedPhase + Physics::pi / 2 + Physics::two_pi, Physics::two_pi);

        if (!isDCGun) {
            INFOMSG(level1 << ">>>>>> APVETO >>>>>> " << endl);
        }
        INFOMSG(level1 << std::fixed << std::setprecision(4)
                << itsCavity_m->getName() << "_phi = "  << originalPhase * Physics::rad2deg <<  " [deg], "
                << "corresp. in Astra = " << AstraPhase * Physics::rad2deg << " [deg],\n"
                << "E = " << finalEnergy << " [MeV], " << "phi_nom = " << originalPhase * Physics::rad2deg << " [deg]\n"
                << "Ez_0 = " << amplitude << " [MV/m]" << "\n"
                << "time = " << (t + tErr) * 1e9 << " [ns], dt = " << dt * 1e12 << " [ps]" << endl);
        if (!isDCGun) {
            INFOMSG(level1 << " <<<<<< APVETO <<<<<< " << endl);
        }

        optimizedPhase = originalPhase;
    }
    INFOMSG(level1 << "* " << std::right << std::setw(83) << std::setfill('*') << "*\n"
            << std::setfill(' ') << endl);

    return optimizedPhase;
}

double CavityAutophaser::guessCavityPhase(double t) {
    const Vector_t &refP = initialP_m;
    double Phimax = 0.0;
    bool apVeto;
    RFCavity *element = static_cast<RFCavity *>(itsCavity_m.get());
    double orig_phi = element->getPhasem();
    apVeto = element->getAutophaseVeto();
    if (apVeto) {
        return orig_phi;
    }

    Phimax = element->getAutoPhaseEstimate(getEnergyMeV(refP),
                                           t,
                                           itsReference_m.getQ(),
                                           itsReference_m.getM() * 1e-6);

    return std::fmod(Phimax + Physics::two_pi, Physics::two_pi);
}

std::pair<double, double> CavityAutophaser::optimizeCavityPhase(double initialPhase,
                                                                double t,
                                                                double dt) {

    RFCavity *element = static_cast<RFCavity *>(itsCavity_m.get());
    double originalPhase = element->getPhasem();

    if (element->getAutophaseVeto()) {
        double basePhase = std::fmod(element->getFrequencym() * t, Physics::two_pi);
        double phase = std::fmod(originalPhase - basePhase + Physics::two_pi, Physics::two_pi);
        double E = track(initialR_m, initialP_m, t, dt, phase);
        std::pair<double, double> status(originalPhase, E);//-basePhase, E);
        return status;
    }

    double Phimax = initialPhase;
    double phi = initialPhase;
    double dphi = Physics::pi / 360.0;
    const int numRefinements = Options::autoPhase;

    int j = -1;
    double E = track(initialR_m, initialP_m, t, dt, phi);
    double Emax = E;

    do {
        j ++;
        Emax = E;
        initialPhase = phi;
        phi -= dphi;
        E = track(initialR_m, initialP_m, t, dt, phi);
    } while(E > Emax);

    if(j == 0) {
        phi = initialPhase;
        E = Emax;
        // j = -1;
        do {
            // j ++;
            Emax = E;
            initialPhase = phi;
            phi += dphi;
            E = track(initialR_m, initialP_m, t, dt, phi);
        } while(E > Emax);
    }

    for(int rl = 0; rl < numRefinements; ++ rl) {
        dphi /= 2.;
        phi = initialPhase - dphi;
        E = track(initialR_m, initialP_m, t, dt, phi);
        if(E > Emax) {
            initialPhase = phi;
            Emax = E;
        } else {
            phi = initialPhase + dphi;
            E = track(initialR_m, initialP_m, t, dt, phi);
            if(E > Emax) {
                initialPhase = phi;
                Emax = E;
            }
        }
    }
    Phimax = std::fmod(initialPhase + Physics::two_pi, Physics::two_pi);

    E = track(initialR_m, initialP_m, t, dt, Phimax + originalPhase);
    std::pair<double, double> status(Phimax, E);

    return status;
}

double CavityAutophaser::track(Vector_t R,
                               Vector_t P,
                               double t,
                               const double dt,
                               const double phase,
                               std::ofstream *out) const {
    const Vector_t &refP = initialP_m;

    RFCavity *rfc = static_cast<RFCavity *>(itsCavity_m.get());
    double initialPhase = rfc->getPhasem();
    rfc->setPhasem(phase);

    std::pair<double, double> pe = rfc->trackOnAxisParticle(refP(2),
                                                            t,
                                                            dt,
                                                            itsReference_m.getQ(),
                                                            itsReference_m.getM() * 1e-6,
                                                            out);
    rfc->setPhasem(initialPhase);

    double finalKineticEnergy = Util::getEnergy(Vector_t(0.0, 0.0, pe.first), itsReference_m.getM() * 1e-6);

    return finalKineticEnergy;
}