//
// Class OrbitThreader
//
// This class determines the design path by tracking the reference particle through
// the 3D lattice.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2022 Christof Metzger-Kraus
//
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

#include "Algorithms/OrbitThreader.h"

#include "AbsBeamline/BendBase.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "BasicActions/Option.h"
#include "BeamlineCore/MarkerRep.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>

#include <cmath>
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
                             StepSizeConfig stepSizes,
                             OpalBeamline &bl) :
    r_m(r),
    p_m(p),
    pathLength_m(s),
    time_m(t),
    dt_m(dt),
    stepSizes_m(stepSizes),
    zstop_m(stepSizes.getFinalZStop() + std::copysign(1.0, dt) * 2 * maxDiffZBunch),
    itsOpalBeamline_m(bl),
    errorFlag_m(0),
    integrator_m(ref),
    reference_m(ref)
{
    auto opal = OpalData::getInstance();
    if (Ippl::myNode() == 0 && !opal->isOptimizerRun()) {
        std::string fileName = Util::combineFilePath({
            opal->getAuxiliaryOutputDirectory(),
            opal->getInputBasename() + "_DesignPath.dat"
        });
        if (opal->getOpenMode() == OpalData::OpenMode::WRITE ||
            !boost::filesystem::exists(fileName)) {
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

        loggingFrequency_m = std::max(1.0, std::round(1e-11 / std::abs(dt_m)));
    } else {
        loggingFrequency_m = std::numeric_limits<size_t>::max();
    }
    pathLengthRange_m = stepSizes_m.getPathLengthRange();
    pathLengthRange_m.enlargeIfOutside(pathLength_m);
    pathLengthRange_m.enlargeIfOutside(zstop_m);
    stepRange_m.enlargeIfOutside(0);
    stepRange_m.enlargeIfOutside(stepSizes_m.getNumStepsFinestResolution());
    distTrackBack_m = std::min(pathLength_m, std::max(0.0, maxDiffZBunch));
    computeBoundingBox();
}

void OrbitThreader::checkElementLengths(const std::set<std::shared_ptr<Component>>& fields) {
    while (!stepSizes_m.reachedEnd() && pathLength_m > stepSizes_m.getZStop()) {
        ++ stepSizes_m;
    }
    if (stepSizes_m.reachedEnd()) {
        return;
    }
    double driftLength = Physics::c * std::abs(stepSizes_m.getdT()) * euclidean_norm(p_m) / Util::getGamma(p_m);
    for (const std::shared_ptr<Component>& field : fields) {
        double length = field->getElementLength();
        int numSteps = field->getRequiredNumberOfTimeSteps();

        if (length < numSteps * driftLength) {
            throw OpalException("OrbitThreader::checkElementLengths",
                                "The time step is too long compared to the length of the\n"
                                "element '" + field->getName() + "'\n" +
                                "The length of the element is: " + std::to_string(length) + "\n"
                                "The distance the particles drift in " + std::to_string(numSteps) +
                                " time step(s) is: " + std::to_string(numSteps * driftLength));
        }
    }
}

void OrbitThreader::execute() {
    double initialPathLength = pathLength_m;

    auto allElements = itsOpalBeamline_m.getElementByType(ElementType::ANY);
    std::set<std::string> visitedElements;

    trackBack();
    updateBoundingBoxWithCurrentPosition();
    pathLengthRange_m.enlargeIfOutside(pathLength_m);

    Vector_t nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR *= Physics::c * dt_m;

    setDesignEnergy(allElements, visitedElements);

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    std::set<std::shared_ptr<Component>> intersection, currentSet;
    errorFlag_m = EVERYTHINGFINE;
    do {
        checkElementLengths(elementSet);
        if (containsCavity(elementSet)) {
            autophaseCavities(elementSet, visitedElements);
        }

        double initialS = pathLength_m;
        Vector_t initialR = r_m;
        Vector_t initialP = p_m;
        double maxDistance = computeDriftLengthToBoundingBox(elementSet, r_m, p_m);
        integrate(elementSet, maxDistance);

        registerElement(elementSet, initialS,  initialR, initialP);

        if (errorFlag_m == HITMATERIAL) {
            // Shouldn't be reached because reference particle
            // isn't stopped by collimators
            pathLength_m += std::copysign(1.0, dt_m);
        }

        imap_m.add(initialS, pathLength_m, elementSet);

        IndexMap::value_t::const_iterator it = elementSet.begin();
        const IndexMap::value_t::const_iterator end = elementSet.end();
        for (; it != end; ++ it) {
            visitedElements.insert((*it)->getName());
        }

        setDesignEnergy(allElements, visitedElements);

        currentSet = elementSet;
        if (errorFlag_m == EVERYTHINGFINE) {
            nextR = r_m / (Physics::c * dt_m);
            integrator_m.push(nextR, p_m, dt_m);
            nextR *= Physics::c * dt_m;

            elementSet = itsOpalBeamline_m.getElements(nextR);
        }
        intersection.clear();
        std::set_intersection(currentSet.begin(), currentSet.end(),
                              elementSet.begin(), elementSet.end(),
                              std::inserter(intersection, intersection.begin()));
    } while (errorFlag_m != HITMATERIAL &&
             errorFlag_m != EOL &&
             stepRange_m.isInside(currentStep_m) &&
             !(pathLengthRange_m.isOutside(pathLength_m) &&
               intersection.empty() && !(elementSet.empty() || currentSet.empty())));

    imap_m.tidyUp(zstop_m);
    *gmsg << level1 << "\n" << imap_m << endl;
    imap_m.saveSDDS(initialPathLength);

    processElementRegister();
}

void OrbitThreader::integrate(const IndexMap::value_t &activeSet, double /*maxDrift*/) {
    CoordinateSystemTrafo labToBeamline = itsOpalBeamline_m.getCSTrafoLab2Local();
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

        if (((pathLength_m > 0.0 &&
              pathLength_m < zstop_m) || dt_m < 0.0) &&
            currentStep_m % loggingFrequency_m == 0 && Ippl::myNode() == 0 &&
            !OpalData::getInstance()->isOptimizerRun()) {
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
                     << std::setw(18) << std::setprecision(8) << reference_m.getM() * (sqrt(dot(p_m, p_m) + 1) - 1) * Units::eV2MeV
                     << std::setw(18) << std::setprecision(8) << (time_m + 0.5 * dt_m) * Units::s2ns
                     << names
                     << std::endl;
        }

        r_m /= Physics::c * dt_m;
        integrator_m.kick(r_m, p_m, Ef, Bf, dt_m);
        integrator_m.push(r_m, p_m, dt_m);
        r_m *= Physics::c * dt_m;

        pathLength_m += std::copysign(euclidean_norm(r_m - oldR), dt_m);
        ++ currentStep_m;
        time_m += dt_m;

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR *= Physics::c * dt_m;

        if (activeSet.empty() &&
            (pathLengthRange_m.isOutside(pathLength_m) ||
             stepRange_m.isOutside(currentStep_m))) {
            errorFlag_m = EOL;
            globalBoundingBox_m.enlargeToContainPosition(r_m);
            return;
        }

    } while (activeSet == itsOpalBeamline_m.getElements(nextR));
}

bool OrbitThreader::containsCavity(const IndexMap::value_t &activeSet) {
    IndexMap::value_t::const_iterator it = activeSet.begin();
    const IndexMap::value_t::const_iterator end = activeSet.end();

    for (; it != end; ++ it) {
        if ((*it)->getType() == ElementType::TRAVELINGWAVE ||
            (*it)->getType() == ElementType::RFCAVITY) {
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
        if (((*it)->getType() == ElementType::TRAVELINGWAVE ||
             (*it)->getType() == ElementType::RFCAVITY) &&
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
        if ((*it)->getType() == ElementType::TRAVELINGWAVE ||
            (*it)->getType() == ElementType::RFCAVITY) {
            const RFCavity *element = static_cast<const RFCavity *>((*it).get());
            designEnergy = std::max(designEnergy, element->getDesignEnergy());
        }
    }

    return designEnergy;
}

void OrbitThreader::trackBack() {
    dt_m *= -1;
    ValueRange<double> tmpRange;
    std::swap(tmpRange, pathLengthRange_m);
    double initialPathLength = pathLength_m;

    Vector_t nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR *= Physics::c * dt_m;

    while (std::abs(initialPathLength - pathLength_m) < distTrackBack_m) {
        auto elementSet = itsOpalBeamline_m.getElements(nextR);

        double maxDrift = computeDriftLengthToBoundingBox(elementSet, r_m, -p_m);
        maxDrift = std::min(maxDrift, distTrackBack_m);
        integrate(elementSet, maxDrift);

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR *= Physics::c * dt_m;
    }
    std::swap(tmpRange, pathLengthRange_m);
    currentStep_m *= -1;
    dt_m *= -1;

    stepRange_m.enlargeIfOutside(currentStep_m);
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

    auto allElements = itsOpalBeamline_m.getElementByType(ElementType::ANY);
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
            !(element->getType() == ElementType::RFCAVITY ||
              element->getType() == ElementType::TRAVELINGWAVE)) {

            element->setDesignEnergy(kineticEnergyeV);
        }
    }
}

void OrbitThreader::computeBoundingBox() {
    FieldList allElements = itsOpalBeamline_m.getElementByType(ElementType::ANY);
    FieldList::iterator it = allElements.begin();
    const FieldList::iterator end = allElements.end();

    for (; it != end; ++ it) {
        if (it->getElement()->getType() == ElementType::MARKER) {
            continue;
        }
        BoundingBox other = it->getBoundingBoxInLabCoords();
        globalBoundingBox_m.enlargeToContainBoundingBox(other);
    }

    updateBoundingBoxWithCurrentPosition();
}


void OrbitThreader::updateBoundingBoxWithCurrentPosition() {
    Vector_t dR = Physics::c * dt_m * p_m / Util::getGamma(p_m);
    for (const Vector_t& pos : {r_m - 10 * dR, r_m + 10 * dR}) {
        globalBoundingBox_m.enlargeToContainPosition(pos);
    }
}

double OrbitThreader::computeDriftLengthToBoundingBox(const std::set<std::shared_ptr<Component>> & elements,
                                                      const Vector_t & position,
                                                      const Vector_t & direction) const {
    if (elements.empty() ||
        (elements.size() == 1 && (*elements.begin())->getType() == ElementType::DRIFT)) {
        boost::optional<Vector_t> intersectionPoint = globalBoundingBox_m.getIntersectionPoint(position, direction);

        return intersectionPoint ? euclidean_norm(intersectionPoint.get() - position): 10.0;
    }

    return std::numeric_limits<double>::max();
}