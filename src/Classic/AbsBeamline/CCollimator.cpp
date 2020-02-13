#include "AbsBeamline/CCollimator.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Structure/LossDataSink.h"

#include <cmath>
#include <fstream>

extern Inform *gmsg;

CCollimator::CCollimator():CCollimator("")
{}

CCollimator::CCollimator(const std::string &name):
    PluginElement(name) {
    setDimensions(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    setGeom(0.0);
}

CCollimator::CCollimator(const CCollimator &right):
    PluginElement(right),
    informed_m(right.informed_m) {
    setDimensions(right.xstart_m, right.xend_m,
                  right.ystart_m, right.yend_m,
                  right.zstart_m, right.zend_m,
                  right.width_m);
    setGeom(width_m);
}

CCollimator::~CCollimator() {}

void CCollimator::accept(BeamlineVisitor &visitor) const {
    visitor.visitCCollimator(*this);
}

bool CCollimator::doPreCheck(PartBunchBase<double, 3> *bunch) {
    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);

    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        // interested in absolute minimum and maximum
        double xmin = std::min(std::abs(rmin(0)), std::abs(rmax(0)));
        double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
        double ymin = std::min(std::abs(rmin(1)), std::abs(rmax(1)));
        double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
        double rbunch_min = std::hypot(xmin, ymin);
        double rbunch_max = std::hypot(xmax, ymax);

        if( rbunch_max > rmin_m && rbunch_min < rmax_m ){ // check similar to z
            return true;
        }
    }
    return false;
}

// rectangle collimators in cyclotron cylindrical coordinates
// when there is no particlematterinteraction, the particle hitting collimator is deleted directly
bool CCollimator::doCheck(PartBunchBase<double, 3> *bunch, const int /*turnnumber*/, const double /*t*/, const double /*tstep*/) {

    bool flagNeedUpdate = false;
    size_t tempnum = bunch->getLocalNum();
    int pflag = 0;
    // now check each particle in bunch
    for (unsigned int i = 0; i < tempnum; ++i) {
        if (bunch->PType[i] == ParticleType::REGULAR && bunch->R[i](2) < zend_m && bunch->R[i](2) > zstart_m ) {
            // only now careful check in r
            pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
            /// bunch->Bin[i] != -1 makes sure the particle is not stored in more than one collimator
            if ((pflag != 0) && (bunch->Bin[i] != -1)) {
                if (!parmatint_m)
                    lossDs_m->addParticle(bunch->R[i], bunch->P[i], bunch->ID[i]);
                bunch->Bin[i] = -1;
                flagNeedUpdate = true;
            }
        }
    }
    return flagNeedUpdate;
}

bool CCollimator::doFinaliseCheck(PartBunchBase<double, 3> *bunch, bool flagNeedUpdate) {
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    if (flagNeedUpdate && parmatint_m) {
        Vector_t rmin, rmax;
        bunch->get_bounds(rmin, rmax);
        std::pair<Vector_t, double> boundingSphere;
        boundingSphere.first = 0.5 * (rmax + rmin);
        boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);
        parmatint_m->apply(bunch, boundingSphere);
    }
    return flagNeedUpdate;
}

void CCollimator::doInitialise(PartBunchBase<double, 3> */*bunch*/) {
    parmatint_m = getParticleMatterInteraction();
}

void CCollimator::goOnline(const double &) {
    print();
    // PosX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // PosY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // PosZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // MomentumX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // MomentumY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // MomentumZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // time_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // id_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    online_m = true;
}

void CCollimator::doFinalise() {
    *gmsg << "* Finalize cyclotron collimator " << getName() << endl;
}

void CCollimator::print() {
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

    *gmsg << "* CCollimator angle start " << xstart_m << " (Deg) angle end " << ystart_m << " (Deg) "
          << "R start " << xend_m << " (mm) R rend " << yend_m << " (mm)" << endl;
}

void CCollimator::setDimensions(double xstart, double xend, double ystart, double yend, double zstart, double zend, double width) {
    setDimensions(xstart, xend, ystart, yend);
    zstart_m = zstart;
    zend_m   = zend;
    width_m  = width;
    // zstart and zend are independent from x, y
    if (zstart_m > zend_m){
        std::swap(zstart_m, zend_m);
    }
    setGeom(width_m);
}

void CCollimator::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementBase::ElementType CCollimator::getType() const {
    return CCOLLIMATOR;
}

void CCollimator::doSetGeom() {
    // calculate maximum r from these
    rmax_m = -std::numeric_limits<double>::max();
    for (int i=0; i<4; i++) {
        double r = std::hypot(geom_m[i].x, geom_m[i].y);
        rmax_m = std::max(rmax_m, r);
    }
    // some debug printout
    // for (int i=0; i<4; i++) {
    //     *gmsg << "point " << i << " ( " << geom_m[i].x << ", " << geom_m[i].y << ")" << endl;
    // }
    // *gmsg << "rmin " << rmin_m << " rmax " << rmax_m << endl;
}

