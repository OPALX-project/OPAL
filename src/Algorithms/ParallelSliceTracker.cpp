#include "Algorithms/ParallelSliceTracker.h"
#include "Algorithms/OrbitThreader.h"
#include "Algorithms/IndexMap.h"
#include "Algorithms/CavityAutophaser.h"

#include "AbstractObjects/OpalData.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Distribution/Distribution.h"
#include "Lines/Sequence.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "AbsBeamline/SpecificElementVisitor.h"

#include <memory>
#include <string>

class PartData;

ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        const PartData &reference, bool revBeam, bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack) {

    CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
                                   beamline.getInitialDirection());
    referenceToLabCSTrafo_m = labToRef.inverted();

}


ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
                                           EnvelopeBunch &bunch,
                                           DataSink &ds,
                                           const PartData &reference,
                                           bool revBeam,
                                           bool revTrack,
                                           const std::vector<unsigned long long> &maxSteps,
                                           double zstart,
                                           const std::vector<double> &zstop,
                                           const std::vector<double> &dt):
    Tracker(beamline, reference, revBeam, revTrack),
    itsOpalBeamline_m(),
    zstart_m(zstart) {

    itsBunch_m      = &bunch;
    itsDataSink_m   = &ds;

    CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
                                   beamline.getInitialDirection());
    referenceToLabCSTrafo_m = labToRef.inverted();

    for (std::vector<unsigned long long>::const_iterator it = maxSteps.begin(); it != maxSteps.end(); ++ it) {
        localTrackSteps_m.push(*it);
    }
    for (std::vector<double>::const_iterator it = dt.begin(); it != dt.end(); ++ it) {
        dtAllTracks_m.push(*it);
    }
    for (std::vector<double>::const_iterator it = zstop.begin(); it != zstop.end(); ++ it) {
        zStop_m.push(*it);
    }

    timeIntegrationTimer1_m  = IpplTimings::getTimer("Time integration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("Time integration2");
    timeFieldEvaluation_m    = IpplTimings::getTimer("Field evaluation");
    BinRepartTimer_m         = IpplTimings::getTimer("Time of Binary repart.");
    WakeFieldTimer_m         = IpplTimings::getTimer("Time of Wake Field calc.");

    cavities_m = itsOpalBeamline_m.getElementByType(ElementBase::RFCAVITY);
    auto travelingwaves = itsOpalBeamline_m.getElementByType(ElementBase::TRAVELINGWAVE);
    cavities_m.merge(travelingwaves, ClassicField::SortAsc);
}


ParallelSliceTracker::~ParallelSliceTracker()
{ }


void ParallelSliceTracker::visitBeamline(const Beamline &bl) { // borrowed from ParallelTTracker
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    if (fbl->getRelativeFlag()) {
        OpalBeamline stash(fbl->getOrigin3D(), fbl->getInitialDirection());
        stash.swap(itsOpalBeamline_m);
        fbl->iterate(*this, false);
        itsOpalBeamline_m.prepareSections();
        itsOpalBeamline_m.compute3DLattice();
        stash.merge(itsOpalBeamline_m);
        stash.swap(itsOpalBeamline_m);
    } else {
        fbl->iterate(*this, false);
    }
}


void ParallelSliceTracker::saveCavityPhases() { // borrowed from ParallelTTracker
    itsDataSink_m->storeCavityInformation();
}

void ParallelSliceTracker::restoreCavityPhases() { // borrowed from ParallelTTracker
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if (OpalData::getInstance()->hasPriorTrack() ||
        OpalData::getInstance()->inRestartRun()) {
        iterator_t it = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for (; it < end; ++ it) {
            updateRFElement((*it).first, (*it).second);
        }
    }
}

/*
 * The maximum phase is added to the nominal phase of
 * the element. This is done on all nodes except node 0 where
 * the Autophase took place.
 */
void ParallelSliceTracker::updateRFElement(std::string elName, double maxPhase) {
    double phase = 0.0;

    for (FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++fit) {
        if ((*fit).getElement()->getName() == elName) {
            if ((*fit).getElement()->getType() == ElementBase::TRAVELINGWAVE) {
                phase  =  static_cast<TravelingWave *>((*fit).getElement().get())->getPhasem();
                static_cast<TravelingWave *>((*fit).getElement().get())->setPhasem(phase + maxPhase);
            } else {
                phase  = static_cast<RFCavity *>((*fit).getElement().get())->getPhasem();
                static_cast<RFCavity *>((*fit).getElement().get())->setPhasem(phase + maxPhase);
            }

            break;
        }
    }
}


/*
 * All RF-Elements gets updated, where the phiShift is the
 * global phase shift in units of seconds.
 */
void ParallelSliceTracker::printRFPhases() {
    Inform msg("ParallelSliceTracker ");

    FieldList &cl = cavities_m;

    const double RADDEG = 180.0 / Physics::pi;
    const double globalTimeShift = OpalData::getInstance()->getGlobalPhaseShift();

    msg << "\n-------------------------------------------------------------------------------------\n";

    for (FieldList::iterator it = cl.begin(); it != cl.end(); ++it) {
        std::shared_ptr<Component> element(it->getElement());
        std::string name = element->getName();
        double frequency;
        double phase;

        if (element->getType() == ElementBase::TRAVELINGWAVE) {
            phase = static_cast<TravelingWave *>(element.get())->getPhasem();
            frequency = static_cast<TravelingWave *>(element.get())->getFrequencym();
        } else {
            phase = static_cast<RFCavity *>(element.get())->getPhasem();
            frequency = static_cast<RFCavity *>(element.get())->getFrequencym();
        }

        msg << (it == cl.begin()? "": "\n")
            << name
            << ": phi = phi_nom + phi_maxE + global phase shift = " << (phase - globalTimeShift * frequency) * RADDEG << " degree, "
            << "(global phase shift = " << -globalTimeShift *frequency *RADDEG << " degree) \n";
    }

    msg << "-------------------------------------------------------------------------------------\n"
        << endl;
}

void ParallelSliceTracker::execute() {

    Inform msg("ParallelSliceTracker", *gmsg);
    BorisPusher pusher(itsReference);
    const double globalTimeShift = itsBunch_m->doEmission()? OpalData::getInstance()->getGlobalPhaseShift(): 0.0;

    currentSimulationTime_m = itsBunch_m->getT();

    setTime();
    setLastStep();

    prepareSections();
    handleAutoPhasing();

    std::queue<double> timeStepSizes(dtAllTracks_m);
    std::queue<unsigned long long> numSteps(localTrackSteps_m);
    double minTimeStep = timeStepSizes.front();
    unsigned long long totalNumSteps = 0;
    while (timeStepSizes.size() > 0) {
        if (minTimeStep > timeStepSizes.front()) {
            totalNumSteps = std::ceil(totalNumSteps * minTimeStep / timeStepSizes.front());
            minTimeStep = timeStepSizes.front();
        }
        totalNumSteps += std::ceil(numSteps.front() * timeStepSizes.front() / minTimeStep);

        numSteps.pop();
        timeStepSizes.pop();
    }

    itsOpalBeamline_m.activateElements();

    if ( OpalData::getInstance()->hasPriorTrack() ||
         OpalData::getInstance()->inRestartRun()) {

        referenceToLabCSTrafo_m = itsBunch_m->toLabTrafo_m;
        RefPartR_m = referenceToLabCSTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
        RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);

        pathLength_m = itsBunch_m->get_sPos();
        zstart_m = pathLength_m;

        restoreCavityPhases();
    } else {
        RefPartR_m = Vector_t(0.0);
        RefPartP_m = euclidean_norm(itsBunch_m->get_pmean_Distribution()) * Vector_t(0, 0, 1);

        if (itsBunch_m->getTotalNum() > 0) {
            if (!itsOpalBeamline_m.containsSource()) {
                RefPartP_m = itsReference.getP() / itsBunch_m->getM() * Vector_t(0, 0, 1);
            }

            if (zstart_m > pathLength_m) {
                findStartPosition(pusher);
            }

            itsBunch_m->set_sPos(pathLength_m);
        }
    }

    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);
    OrbitThreader oth(itsReference,
                      referenceToLabCSTrafo_m.transformTo(RefPartR_m),
                      referenceToLabCSTrafo_m.rotateTo(RefPartP_m),
                      pathLength_m,
                      0.0,//-rmin(2),
                      itsBunch_m->getT(),
                      minTimeStep,
                      totalNumSteps,
                      zStop_m.back() + 2 * rmax(2),
                      itsOpalBeamline_m);

    oth.execute();

    saveCavityPhases();

    setTime();

    double t = itsBunch_m->getT() - globalTimeShift;
    itsBunch_m->setT(t);

    itsBunch_m->RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
    itsBunch_m->RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);

    *gmsg << level1 << *itsBunch_m << endl;

    unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << "Track start at: " << myt1.time() << ", t= " << Util::getTimeString(t) << "; "
          << "zstart at: " << Util::getLengthString(pathLength_m)
          << endl;

    // prepareEmission();

    *gmsg << level1
          << "Executing ParallelSliceTracker, initial dt= " << Util::getTimeString(itsBunch_m->getdT()) << ";\n"
          << "max integration steps " << getMaxSteps(localTrackSteps_m) << ", next step= " << step << endl
          << "the mass is: " << itsReference.getM() * 1e-6 << " MeV, its charge: "
          << itsReference.getQ() << " [e]" << endl;


    // setOptionalVariables();
    *gmsg << std::scientific;
    itsBunch_m->toLabTrafo_m = referenceToLabCSTrafo_m;
    while (localTrackSteps_m.size() > 0) {
        localTrackSteps_m.front() += step;
        dtCurrentTrack_m = dtAllTracks_m.front();
        selectDT();

        for (; step < localTrackSteps_m.front(); ++step) {

            globalEOL_m = false;

            switchElements();
            computeExternalFields(oth);

            //reduce(&globalEOL_m, &globalEOL_m, OpBitwiseOrAssign());
            //reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
            MPI_Allreduce(MPI_IN_PLACE, &globalEOL_m, 1, MPI_INT, MPI_LAND, Ippl::getComm());

            computeSpaceChargeFields();
            timeIntegration();

            if (step % 1000 + 1 == 1000) {
                msg << level1;
            } else if (step % 100 + 1 == 100) {
                msg << level2;
            } else {
                msg << level3;
            }

            msg << " Step " << step << " at " << itsBunch_m->zAvg() << " [m] t= "
                << itsBunch_m->getT() << " [s] E=" << itsBunch_m->Eavg() * 1e-6
                << " [MeV]" << endl;

            currentSimulationTime_m += itsBunch_m->getdT();
            itsBunch_m->setT(currentSimulationTime_m);

            dumpStats(step);

            if (hasEndOfLineReached())
                break;

            double gamma = itsBunch_m->Eavg() / (Physics::m_e * 1e9) + 1.0;
            double beta = sqrt(1.0 - 1.0 / (gamma * gamma));

            double driftPerTimeStep = itsBunch_m->getdT() * Physics::c * beta;
            if (std::abs(zStop_m.front() - itsBunch_m->zAvg()) < 0.5 * driftPerTimeStep)
                localTrackSteps_m.front() = step;
        }

        if (globalEOL_m)
            break;

        dtAllTracks_m.pop();
        localTrackSteps_m.pop();
        zStop_m.pop();
    }

    OpalData::getInstance()->setLastStep(step);
    writeLastStepPhaseSpace(step, itsBunch_m->get_sPos());
    itsOpalBeamline_m.switchElementsOff();
    msg << "done executing ParallelSliceTracker" << endl;
}


void ParallelSliceTracker::timeIntegration() {

    IpplTimings::startTimer(timeIntegrationTimer1_m);
    itsBunch_m->timeStep(itsBunch_m->getdT());
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}


void ParallelSliceTracker::handleAutoPhasing() {

    if (Options::autoPhase == 0) return;

    if(!OpalData::getInstance()->inRestartRun()) {
        itsDataSink_m->storeCavityInformation();
    }

    auto it = OpalData::getInstance()->getFirstMaxPhases();
    auto end = OpalData::getInstance()->getLastMaxPhases();
    for (; it != end; ++ it) {
        updateRFElement((*it).first, (*it).second);
    }
    printRFPhases();
}


void ParallelSliceTracker::prepareSections() {

    itsBeamline_m.accept(*this);
    itsOpalBeamline_m.prepareSections();

    itsOpalBeamline_m.compute3DLattice();
    itsOpalBeamline_m.save3DLattice();
    itsOpalBeamline_m.save3DInput();
}


void ParallelSliceTracker::computeExternalFields(OrbitThreader &oth) {

    IpplTimings::startTimer(timeFieldEvaluation_m);

    Vector_t externalE, externalEThis, externalB, externalBThis, KR, KRThis, KT, KTThis;

    const unsigned int localNum = itsBunch_m->getLocalNum();
    // bool locPartOutOfBounds = false, globPartOutOfBounds = false;
    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);
    IndexMap::value_t elements;

    try {
        elements = oth.query(pathLength_m + 0.5 * (rmax(2) + rmin(2)), rmax(2) - rmin(2));
    } catch(IndexMap::OutOfBounds &e) {
        globalEOL_m = true;
        return;
    }

    const double time = itsBunch_m->getT() + 0.5 * itsBunch_m->getdT();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (unsigned int i = 0; i < localNum; ++ i) {

        externalB     = Vector_t(0.0);
        externalE     = Vector_t(0.0);
        KR            = Vector_t(0.0);
        KT            = Vector_t(0.0);

        IndexMap::value_t::const_iterator it = elements.begin();
        for (; it != end; ++ it) {
            const CoordinateSystemTrafo toLocal = (itsOpalBeamline_m.getMisalignment((*it)) *
                                                   itsOpalBeamline_m.getCSTrafoLab2Local((*it)));
            const CoordinateSystemTrafo fromLocal = toLocal.inverted();

            //FIXME: why not x=y=0.0?
            // for (unsigned int k = 0; k < 100; ++ k) {
            externalEThis = Vector_t(0.0);
            externalBThis = Vector_t(0.0);
            KRThis        = Vector_t(0.0);
            KTThis        = Vector_t(0.0);

            Vector_t pos = itsBunch_m->getR(i);
            Vector_t mom = itsBunch_m->getP(i);

            itsBunch_m->setR(i, toLocal.transformTo(pos));
            itsBunch_m->setP(i, toLocal.rotateTo(mom));

            if ((*it)->apply(itsBunch_m->getR(i),
                             itsBunch_m->getP(i),
                             time,
                             externalEThis,
                             externalBThis)) {
                itsBunch_m->Bin[i] = -1;
                continue;
            }

            (*it)->addKR(i, time, KRThis);
            (*it)->addKT(i, time, KTThis);

            externalE += fromLocal.rotateTo(externalEThis);
            externalB += fromLocal.rotateTo(externalBThis);

            KR += fromLocal.rotateTo(KRThis);
            KT += fromLocal.rotateTo(KTThis);

            itsBunch_m->setR(i, pos);
            itsBunch_m->setP(i, mom);
        }

        itsBunch_m->setExternalFields(i, externalE, externalB, KR, KT);
    }

    IpplTimings::stopTimer(timeFieldEvaluation_m);
}


void ParallelSliceTracker::computeSpaceChargeFields() {

    itsBunch_m->computeSpaceCharge();
}


void ParallelSliceTracker::dumpStats(long long step) {

    double sposRef = itsBunch_m->get_sPos();
    if (step != 0 && (step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0))
        writePhaseSpace(step, sposRef);
}


void ParallelSliceTracker::switchElements(double /*scaleMargin*/) {

    double margin = 1.0;

    currentSimulationTime_m = itsBunch_m->getT();
    itsOpalBeamline_m.switchElements(itsBunch_m->zTail() - margin,
                                      itsBunch_m->zHead() + margin,
                                      itsBunch_m->Eavg() * 1e-6);
}


void ParallelSliceTracker::setLastStep() {

    // unsigned long long step = 0;

    // if (OpalData::getInstance()->inRestartRun()) {

    //     int prevDumpFreq = OpalData::getInstance()->getRestartDumpFreq();
    //     step = OpalData::getInstance()->getRestartStep() * prevDumpFreq + 1;
    //     maxSteps_m += step;
    // } else {

    //     step = OpalData::getInstance()->getLastStep() + 1;
    //     maxSteps_m += step;
    // }

    // OpalData::getInstance()->setLastStep(step);
}


void ParallelSliceTracker::setTime() {

    if (OpalData::getInstance()->inRestartRun())
        currentSimulationTime_m = itsBunch_m->getT();
    else
        currentSimulationTime_m = itsBunch_m->getT();
}


bool ParallelSliceTracker::hasEndOfLineReached() {
    reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
    return globalEOL_m;
}

void ParallelSliceTracker::findStartPosition(const BorisPusher &pusher) { // borrowed from ParallelTTracker

    double t = 0.0;
    itsBunch_m->setT(t);

    dtCurrentTrack_m = dtAllTracks_m.front();
    selectDT();

    if (Util::getEnergy(RefPartP_m, itsBunch_m->getM()) < 1e-3) {
        double gamma = 0.1 / itsBunch_m->getM() + 1.0;
        RefPartP_m = sqrt(std::pow(gamma, 2) - 1) * Vector_t(0, 0, 1);
    }

    while (true) {
        autophaseCavities(pusher);

        t += itsBunch_m->getdT();
        itsBunch_m->setT(t);

        Vector_t oldR = RefPartR_m;
        updateReferenceParticle(pusher);
        pathLength_m += euclidean_norm(RefPartR_m - oldR);

        if (pathLength_m > zStop_m.front()) {
            if (localTrackSteps_m.size() == 0) return;

            dtAllTracks_m.pop();
            localTrackSteps_m.pop();
            zStop_m.pop();

            selectDT();
        }

        double speed = euclidean_norm(RefPartP_m) * Physics::c / sqrt(dot(RefPartP_m, RefPartP_m) + 1);
        if (std::abs(pathLength_m - zstart_m) <=  0.5 * itsBunch_m->getdT() * speed) {
            double tau = (pathLength_m - zstart_m) / speed;

            t += tau;
            itsBunch_m->setT(t);

            RefPartR_m /= (Physics::c * tau);
            pusher.push(RefPartR_m, RefPartP_m, tau);
            RefPartR_m *= (Physics::c * tau);

            pathLength_m = zstart_m;

            CoordinateSystemTrafo update(RefPartR_m,
                                         getQuaternion(RefPartP_m, Vector_t(0, 0, 1)));
            referenceToLabCSTrafo_m = referenceToLabCSTrafo_m * update.inverted();

            RefPartR_m = update.transformTo(RefPartR_m);
            RefPartP_m = update.rotateTo(RefPartP_m);

            return;
        }
    }
}

void ParallelSliceTracker::updateReferenceParticle(const BorisPusher &pusher) { // borrowed from ParallelTTracker
    //static size_t step = 0;
    const double dt = std::min(itsBunch_m->getT(), itsBunch_m->getdT());
    const double scaleFactor = Physics::c * dt;
    Vector_t Ef(0.0), Bf(0.0);
    // Vector_t oldR = RefPartR_m;

    RefPartR_m /= scaleFactor;
    pusher.push(RefPartR_m, RefPartP_m, dt);
    RefPartR_m *= scaleFactor;

    IndexMap::value_t elements = itsOpalBeamline_m.getElements(referenceToLabCSTrafo_m.transformTo(RefPartR_m));
    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        CoordinateSystemTrafo refToLocalCSTrafo = itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * referenceToLabCSTrafo_m;

        Vector_t localR = refToLocalCSTrafo.transformTo(RefPartR_m);
        Vector_t localP = refToLocalCSTrafo.rotateTo(RefPartP_m);
        Vector_t localE(0.0), localB(0.0);

        if ((*it)->applyToReferenceParticle(localR,
                                            localP,
                                            itsBunch_m->getT() - 0.5 * dt,
                                            localE,
                                            localB)) {
            *gmsg << level1 << "The reference particle hit an element" << endl;
            globalEOL_m = true;
        }

        Ef += refToLocalCSTrafo.rotateFrom(localE);
        Bf += refToLocalCSTrafo.rotateFrom(localB);
    }

    pusher.kick(RefPartR_m, RefPartP_m, Ef, Bf, dt);

    RefPartR_m /= scaleFactor;
    pusher.push(RefPartR_m, RefPartP_m, dt);
    RefPartR_m *= scaleFactor;
    //++ step;
}

void ParallelSliceTracker::autophaseCavities(const BorisPusher &pusher) { // borrowed from ParallelTTracker

    double t = itsBunch_m->getT();
    Vector_t nextR = RefPartR_m / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, RefPartP_m, itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(referenceToLabCSTrafo_m.transformTo(nextR));
    for (auto element: elementSet) {
        if (element->getType() == ElementBase::TRAVELINGWAVE) {
            const TravelingWave *TWelement = static_cast<const TravelingWave *>(element.get());
            if (!TWelement->getAutophaseVeto()) {
                RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, RefPartP_m),
                                       t, itsBunch_m->getdT());
                RefPartR_m = referenceToLabCSTrafo_m.transformFrom(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(RefPartP_m);
            }

        } else if (element->getType() == ElementBase::RFCAVITY) {
            const RFCavity *RFelement = static_cast<const RFCavity *>(element.get());
            if (!RFelement->getAutophaseVeto()) {
                RefPartR_m = referenceToLabCSTrafo_m.transformTo(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateTo(RefPartP_m);
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, RefPartP_m),
                                       t, itsBunch_m->getdT());
                RefPartR_m = referenceToLabCSTrafo_m.transformFrom(RefPartR_m);
                RefPartP_m = referenceToLabCSTrafo_m.rotateFrom(RefPartP_m);
            }
        }
    }
}

void ParallelSliceTracker::selectDT() { // borrowed from ParallelTTracker

    if (itsBunch_m->getIfBeamEmitting()) {
        double dt = itsBunch_m->getEmissionDeltaT();
        itsBunch_m->setdT(dt);
    } else {
        double dt = dtCurrentTrack_m;
        itsBunch_m->setdT(dt);
    }
}

void ParallelSliceTracker::changeDT() { // borrowed from ParallelTTracker
    selectDT();
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++ i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}