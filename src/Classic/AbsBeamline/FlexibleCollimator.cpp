//
// Class FlexibleCollimator
//   Defines the abstract interface for a collimator.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/FlexibleCollimator.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Solvers/ParticleMatterInteractionHandler.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <memory>

extern Inform *gmsg;

FlexibleCollimator::FlexibleCollimator():
    FlexibleCollimator("")
{}


FlexibleCollimator::FlexibleCollimator(const FlexibleCollimator& right):
    Component(right),
    description_m(right.description_m),
    bb_m(right.bb_m),
    tree_m(/*right.tree_m*/),
    informed_m(right.informed_m),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(nullptr)
{
    for (const std::shared_ptr<mslang::Base>& obj: right.holes_m) {
        holes_m.emplace_back(obj->clone());
    }

    tree_m.bb_m = bb_m;
    tree_m.objects_m.insert(tree_m.objects_m.end(), holes_m.begin(), holes_m.end());
    tree_m.buildUp();
}


FlexibleCollimator::FlexibleCollimator(const std::string& name):
    Component(name),
    description_m(""),
    informed_m(false),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(nullptr)
{}


FlexibleCollimator::~FlexibleCollimator() {
    if (online_m)
        goOffline();
    // for (mslang::Base *obj: holes_m) {
    //     delete obj;
    // }
}


void FlexibleCollimator::accept(BeamlineVisitor& visitor) const {
    visitor.visitFlexibleCollimator(*this);
}

bool FlexibleCollimator::isStopped(const Vector_t& R) {
    const double z = R(2);

    if ((z < 0.0) ||
        (z > getElementLength()) ||
        (!isInsideTransverse(R))) {
        return false;
    }

    if (!bb_m.isInside(R)) {
        return getFlagDeleteOnTransverseExit();
    }

    if (!tree_m.isInside(R)) {
        return true;
    }

    return false;
}

bool FlexibleCollimator::apply(const size_t& i, const double& t,
                               Vector_t& /*E*/, Vector_t& /*B*/) {
    const Vector_t& R = RefPartBunch_m->R[i];
    bool pdead = isStopped(R);

    if (pdead) {
        if (lossDs_m) {
            const Vector_t& P = RefPartBunch_m->P[i];
            const double& dt = RefPartBunch_m->dt[i];
            const Vector_t singleStep = Physics::c * dt * Util::getBeta(P);
            double frac = -R(2) / singleStep(2);
            lossDs_m->addParticle(OpalParticle(RefPartBunch_m->ID[i],
                                               R + frac * singleStep, P,
                                               t + frac * dt,
                                               RefPartBunch_m->Q[i], RefPartBunch_m->M[i]));
        }
        ++losses_m;
    }
    return pdead;
}

bool FlexibleCollimator::applyToReferenceParticle(const Vector_t& /*R*/,
                                                  const Vector_t& /*P*/,
                                                  const double& /*t*/,
                                                  Vector_t& /*E*/,
                                                  Vector_t& /*B*/) {
    return false;
}

// rectangle collimators in cyclotron cyclindral coordinates
// without particlematterinteraction, the particle hitting collimator is deleted directly
bool FlexibleCollimator::checkCollimator(PartBunchBase<double, 3>* /*bunch*/,
                                         const int /*turnnumber*/,
                                         const double /*t*/,
                                         const double /*tstep*/) {
    return false;
}

void FlexibleCollimator::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    endField = startField + getElementLength();

    parmatint_m = getParticleMatterInteraction();

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getOutputFN(), !Options::asciidump));

    goOnline(-1e6);
}

void FlexibleCollimator::initialise(PartBunchBase<double, 3>* bunch) {
    RefPartBunch_m = bunch;

    parmatint_m = getParticleMatterInteraction();

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getOutputFN(), !Options::asciidump));

    goOnline(-1e6);
}

void FlexibleCollimator::finalise() {
    if (online_m)
        goOffline();
    *gmsg << "* Finalize flexible collimator " << getName() << endl;
}

void FlexibleCollimator::goOnline(const double&) {
    print();
    online_m = true;
}

void FlexibleCollimator::print() {
    if (RefPartBunch_m == nullptr) {
        if (!informed_m) {
            std::string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
            ERRORMSG(errormsg << endl);
            if (Ippl::myNode() == 0) {
                std::ofstream omsg("errormsg.txt", std::ios_base::app);
                omsg << errormsg << std::endl;
                omsg.close();
            }
            informed_m = true;
        }
        return;
    }

    *gmsg << level3;
}

void FlexibleCollimator::goOffline() {
    if (online_m && lossDs_m)
        lossDs_m->save();
    lossDs_m.reset(0);
    online_m = false;
}

void FlexibleCollimator::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementType FlexibleCollimator::getType() const {
    return ElementType::FLEXIBLECOLLIMATOR;
}

void FlexibleCollimator::setDescription(const std::string& desc) {
    tree_m.reset();
    holes_m.clear();

    mslang::Function* fun;

    if (!mslang::parse(desc, fun))
        throw GeneralClassicException("FlexibleCollimator::setDescription",
                                      "Couldn't parse input file");

    fun->apply(holes_m);

    if (holes_m.size() == 0) return;

    for (std::shared_ptr<mslang::Base>& it: holes_m) {
        it->computeBoundingBox();
    }

    std::shared_ptr<mslang::Base>& first = holes_m.front();
    const mslang::BoundingBox2D& bb = first->bb_m;

    Vector_t llc(bb.center_m[0] - 0.5 * bb.width_m,
                 bb.center_m[1] - 0.5 * bb.height_m,
                 0.0);
    Vector_t urc(bb.center_m[0] + 0.5 * bb.width_m,
                 bb.center_m[1] + 0.5 * bb.height_m,
                 0.0);

    for (const std::shared_ptr<mslang::Base>& it: holes_m) {
        const mslang::BoundingBox2D& bb = it->bb_m;
        llc[0] = std::min(llc[0], bb.center_m[0] - 0.5 * bb.width_m);
        llc[1] = std::min(llc[1], bb.center_m[1] - 0.5 * bb.height_m);
        urc[0] = std::max(urc[0], bb.center_m[0] + 0.5 * bb.width_m);
        urc[1] = std::max(urc[1], bb.center_m[1] + 0.5 * bb.height_m);
    }

    double width = urc[0] - llc[0];
    double height = urc[1] - llc[1];

    llc[0] -= 1e-3 * width;
    urc[0] += 1e-3 * width;
    llc[1] -= 1e-3 * height;
    urc[1] += 1e-3 * height;

    bb_m = mslang::BoundingBox2D(llc, urc);

    tree_m.bb_m = bb_m;
    tree_m.objects_m.insert(tree_m.objects_m.end(), holes_m.begin(), holes_m.end());
    tree_m.buildUp();

    delete fun;
}

void FlexibleCollimator::writeHolesAndQuadtree(const std::string& baseFilename) const {
    if (Ippl::myNode() == 0) {
        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            baseFilename
        });

        std::ofstream out(fname + "_quadtree.gpl");
        tree_m.writeGnuplot(out);
        out.close();

        out.open(fname + "_holes.gpl");
        for (const std::shared_ptr<mslang::Base> &obj: holes_m) {
            obj->writeGnuplot(out);
        }
        out.close();
    }
}
