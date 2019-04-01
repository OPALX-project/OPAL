#include "Algorithms/OrbitThreader.h"
#include "Algorithms/CavityAutophaser.h"

#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"

#include "AbstractObjects/OpalData.h"
#include "BasicActions/Option.h"
#include "Utilities/Options.h"

#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>

#include <iostream>
#include <limits>

#define HITMATERIAL 0x80000000
#define EOL 0x40000000
#define EVERYTHINGFINE 0x00000000
extern Inform *gmsg;

OrbitThreader::OrbitThreader(const PartData &ref,
                             const Vector_t &r,
                             const Vector_t &p,
                             double s,
                             double maxDiffZBunch,
                             double t,
                             double dt,
                             size_t maxIntegSteps,
                             double zstop,
                             OpalBeamline &bl) :
    r_m(r),
    p_m(p),
    pathLength_m(s),
    time_m(t),
    initialTime_m(t),
    dt_m(dt),
    maxIntegSteps_m(maxIntegSteps),
    zstop_m(zstop),
    itsOpalBeamline_m(bl),
    errorFlag_m(0),
    integrator_m(ref),
    reference_m(ref)
{
    auto opal = OpalData::getInstance();
    if (Ippl::myNode() == 0 && !opal->isOptimizerRun()) {
        std::string fileName = "data/" + OpalData::getInstance()->getInputBasename() + "_DesignPath.dat";
        if (Options::openMode == Options::WRITE || !boost::filesystem::exists(fileName)) {
            logger_m.open(fileName);
            logger_m << "#"
                     << std::setw(17) << "1 - s"
                     << std::setw(18) << "2 - Rx"
                     << std::setw(18) << "3 - Ry"
                     << std::setw(18) << "4 - Rz"
                     << std::setw(18) << "5 - Px"
                     << std::setw(18) << "6 - Py"
                     << std::setw(18) << "7 - Pz"
                     << std::setw(18) << "8 - Efx"
                     << std::setw(18) << "9 - Efy"
                     << std::setw(18) << "10 - Efz"
                     << std::setw(18) << "11 - Bfx"
                     << std::setw(18) << "12 - Bfy"
                     << std::setw(18) << "13 - Bfz"
                     << std::setw(18) << "14 - Ekin"
                     << std::setw(18) << "15 - t"
                     << std::endl;
        } else {
            logger_m.open(fileName, std::ios_base::app);
        }

        loggingFrequency_m = std::max(1.0, floor(1e-11 / dt_m + 0.5));
    } else {
        loggingFrequency_m = std::numeric_limits<size_t>::max();
    }

    dZ_m = std::min(pathLength_m, std::max(0.0, maxDiffZBunch));
}

void OrbitThreader::execute() {
    double initialPathLength = pathLength_m;
    trackBack();
    Vector_t nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR *= Physics::c * dt_m;
    auto allElements = itsOpalBeamline_m.getElementByType(ElementBase::ANY);
    std::set<std::string> visitedElements;
    setDesignEnergy(allElements, visitedElements);

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    errorFlag_m = EVERYTHINGFINE;
    do {
        if (containsCavity(elementSet)) {
            autophaseCavities(elementSet, visitedElements);
        }

        double initialS = pathLength_m;
        Vector_t initialR = r_m;
        Vector_t initialP = p_m;
        integrate(elementSet, maxIntegSteps_m);

        registerElement(elementSet, initialS,  initialR, initialP);

        if (errorFlag_m == HITMATERIAL) {
            pathLength_m += 1.0;
        }

        imap_m.add(initialS, pathLength_m, elementSet);

        IndexMap::value_t::const_iterator it = elementSet.begin();
        const IndexMap::value_t::const_iterator end = elementSet.end();
        for (; it != end; ++ it) {
            visitedElements.insert((*it)->getName());
        }

        setDesignEnergy(allElements, visitedElements);

        if (errorFlag_m == EVERYTHINGFINE) {
            nextR = r_m / (Physics::c * dt_m);
            integrator_m.push(nextR, p_m, dt_m);
            nextR *= Physics::c * dt_m;

            elementSet = itsOpalBeamline_m.getElements(nextR);
        }
    } while ((time_m - initialTime_m) / dt_m < maxIntegSteps_m &&
             errorFlag_m != HITMATERIAL &&
             errorFlag_m != EOL);

    // imap_m.tidyUp();
    *gmsg << level1 << "\n" << imap_m << endl;
    imap_m.saveSDDS(initialPathLength);

    processElementRegister();
}

void OrbitThreader::integrate(const IndexMap::value_t &activeSet, size_t maxSteps) {
    static size_t step = 0;
    Vector_t nextR;
    do {
        errorFlag_m = EVERYTHINGFINE;

        IndexMap::value_t::const_iterator it = activeSet.begin();
        const IndexMap::value_t::const_iterator end = activeSet.end();
        Vector_t oldR = r_m;

        r_m /= Physics::c * dt_m;
        integrator_m.push(r_m, p_m, dt_m);
        r_m *= Physics::c * dt_m;

        Vector_t Ef(0.0), Bf(0.0);
        std::string names("\t");
        for (; it != end; ++ it) {
            Vector_t localR = itsOpalBeamline_m.transformToLocalCS(*it, r_m);
            Vector_t localP = itsOpalBeamline_m.rotateToLocalCS(*it, p_m);
            Vector_t localE(0.0), localB(0.0);

            if ((*it)->applyToReferenceParticle(localR, localP, time_m + 0.5 * dt_m, localE, localB)) {
                errorFlag_m = HITMATERIAL;
                return;
            }
            names += (*it)->getName() + ", ";

            Ef += itsOpalBeamline_m.rotateFromLocalCS(*it, localE);
            Bf += itsOpalBeamline_m.rotateFromLocalCS(*it, localB);
        }

        if (step % loggingFrequency_m == 0 && Ippl::myNode() == 0 && !OpalData::getInstance()->isOptimizerRun()) {
            logger_m << std::setw(18) << std::setprecision(8) << pathLength_m + std::copysign(euclidean_norm(r_m - oldR), dt_m)
                     << std::setw(18) << std::setprecision(8) << r_m(0)
                     << std::setw(18) << std::setprecision(8) << r_m(1)
                     << std::setw(18) << std::setprecision(8) << r_m(2)
                     << std::setw(18) << std::setprecision(8) << p_m(0)
                     << std::setw(18) << std::setprecision(8) << p_m(1)
                     << std::setw(18) << std::setprecision(8) << p_m(2)
                     << std::setw(18) << std::setprecision(8) << Ef(0)
                     << std::setw(18) << std::setprecision(8) << Ef(1)
                     << std::setw(18) << std::setprecision(8) << Ef(2)
                     << std::setw(18) << std::setprecision(8) << Bf(0)
                     << std::setw(18) << std::setprecision(8) << Bf(1)
                     << std::setw(18) << std::setprecision(8) << Bf(2)
                     << std::setw(18) << std::setprecision(8) << reference_m.getM() * (sqrt(dot(p_m, p_m) + 1) - 1) * 1e-6
                     << std::setw(18) << std::setprecision(8) << (time_m + 0.5 * dt_m) * 1e9
                     << names
                     << std::endl;
        }

        r_m /= Physics::c * dt_m;
        integrator_m.kick(r_m, p_m, Ef, Bf, dt_m);
        integrator_m.push(r_m, p_m, dt_m);
        r_m *= Physics::c * dt_m;

        pathLength_m += std::copysign(euclidean_norm(r_m - oldR), dt_m);
        ++ step;
        time_m += dt_m;

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR *= Physics::c * dt_m;

        if (pathLength_m > zstop_m) {
            errorFlag_m = EOL;
            return;
        }

    } while (((time_m - initialTime_m) / dt_m < maxSteps) &&
             (activeSet == itsOpalBeamline_m.getElements(nextR)));
}

bool OrbitThreader::containsCavity(const IndexMap::value_t &activeSet) {
    IndexMap::value_t::const_iterator it = activeSet.begin();
    const IndexMap::value_t::const_iterator end = activeSet.end();

    for (; it != end; ++ it) {
        if ((*it)->getType() == ElementBase::TRAVELINGWAVE ||
            (*it)->getType() == ElementBase::RFCAVITY) {
            return true;
        }
    }

    return false;
}

void OrbitThreader::autophaseCavities(const IndexMap::value_t &activeSet,
                                      const std::set<std::string> &visitedElements) {
    if (Options::autoPhase == 0) return;

    IndexMap::value_t::const_iterator it = activeSet.begin();
    const IndexMap::value_t::const_iterator end = activeSet.end();

    for (; it != end; ++ it) {
        if (((*it)->getType() == ElementBase::TRAVELINGWAVE ||
             (*it)->getType() == ElementBase::RFCAVITY) &&
            visitedElements.find((*it)->getName()) == visitedElements.end()) {

            Vector_t initialR = itsOpalBeamline_m.transformToLocalCS(*it, r_m);
            Vector_t initialP = itsOpalBeamline_m.rotateToLocalCS(*it, p_m);

            CavityAutophaser ap(reference_m, *it);
            ap.getPhaseAtMaxEnergy(initialR,
                                   initialP,
                                   time_m,
                                   dt_m);
        }
    }
}


double OrbitThreader::getMaxDesignEnergy(const IndexMap::value_t &elementSet) const
{
    IndexMap::value_t::const_iterator it = elementSet.begin();
    const IndexMap::value_t::const_iterator end = elementSet.end();

    double designEnergy = 0.0;
    for (; it != end; ++ it) {
        if ((*it)->getType() == ElementBase::TRAVELINGWAVE ||
            (*it)->getType() == ElementBase::RFCAVITY) {
            const RFCavity *element = static_cast<const RFCavity *>((*it).get());
            designEnergy = std::max(designEnergy, element->getDesignEnergy());
        }
    }

    return designEnergy;
}

void OrbitThreader::trackBack() {
    dt_m *= -1;
    double initialPathLength = pathLength_m;

    Vector_t nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR *= Physics::c * dt_m;

    while (initialPathLength - pathLength_m < dZ_m) {
        auto elementSet = itsOpalBeamline_m.getElements(nextR);

        integrate(elementSet, 1000);

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR *= Physics::c * dt_m;
    }

    dt_m *= -1;
}

void OrbitThreader::registerElement(const IndexMap::value_t &elementSet,
                                    double start,
                                    const Vector_t &R,
                                    const Vector_t &P) {

    IndexMap::value_t::const_iterator it = elementSet.begin();
    const IndexMap::value_t::const_iterator end = elementSet.end();

    for (; it != end; ++ it) {
        bool found = false;
        std::string name = (*it)->getName();
        auto prior = elementRegistry_m.equal_range(name);

        for (auto pit = prior.first; pit != prior.second; ++ pit) {
            if (std::abs((*pit).second.endField_m - start) < 1e-10) {
                found = true;
                (*pit).second.endField_m = pathLength_m;
                break;
            }
        }

        if (found) continue;

        Vector_t initialR = itsOpalBeamline_m.transformToLocalCS(*it, R);
        Vector_t initialP = itsOpalBeamline_m.rotateToLocalCS(*it, P);
        double elementEdge = start - initialR(2) * euclidean_norm(initialP) / initialP(2);

        elementPosition ep = {start, pathLength_m, elementEdge};
        elementRegistry_m.insert(std::make_pair(name, ep));
    }
}

void OrbitThreader::processElementRegister() {
    std::map<std::string, std::set<elementPosition, elementPositionComp> > tmpRegistry;

    for (auto it = elementRegistry_m.begin(); it != elementRegistry_m.end(); ++ it) {
        const std::string &name = (*it).first;
        elementPosition &ep = (*it).second;

        auto prior = tmpRegistry.find(name);
        if (prior == tmpRegistry.end()) {
            std::set<elementPosition, elementPositionComp> tmpSet;
            tmpSet.insert(ep);
            tmpRegistry.insert(std::make_pair(name, tmpSet));
            continue;
        }

        std::set<elementPosition, elementPositionComp> &set = (*prior).second;
        set.insert(ep);
    }

    auto allElements = itsOpalBeamline_m.getElementByType(ElementBase::ANY);
    FieldList::iterator it = allElements.begin();
    const FieldList::iterator end = allElements.end();
    for (; it != end; ++ it) {
        std::string name = (*it).getElement()->getName();
        auto trit = tmpRegistry.find(name);
        if (trit == tmpRegistry.end()) continue;

        std::queue<std::pair<double, double> > range;
        std::set<elementPosition, elementPositionComp> &set = (*trit).second;

        for (auto sit = set.begin(); sit != set.end(); ++ sit) {
            range.push(std::make_pair((*sit).elementEdge_m, (*sit).endField_m));
        }
        (*it).getElement()->setActionRange(range);
    }
}

void OrbitThreader::setDesignEnergy(FieldList &allElements, const std::set<std::string> &visitedElements) {
    double kineticEnergyeV = reference_m.getM() * (sqrt(dot(p_m, p_m) + 1.0) - 1.0);

    FieldList::iterator it = allElements.begin();
    const FieldList::iterator end = allElements.end();
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        if (visitedElements.find(element->getName()) == visitedElements.end() &&
            !(element->getType() == ElementBase::RFCAVITY ||
              element->getType() == ElementBase::TRAVELINGWAVE)) {

            element->setDesignEnergy(kineticEnergyeV);
        }
    }
}
