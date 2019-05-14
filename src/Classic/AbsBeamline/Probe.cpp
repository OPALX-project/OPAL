#include "AbsBeamline/Probe.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Structure/PeakFinder.h"

extern Inform *gmsg;

Probe::Probe():Probe("")
{}

Probe::Probe(const std::string &name):
    PluginElement(name),
    step_m(0.0)
{}

Probe::Probe(const Probe &right):
    PluginElement(right),
    step_m(right.step_m)
{}

Probe::~Probe() {}

void Probe::accept(BeamlineVisitor &visitor) const {
    visitor.visitProbe(*this);
}

void Probe::doInitialise(PartBunchBase<double, 3> *bunch) {
    *gmsg << "* Initialize probe" << endl;
    bool singlemode = (bunch->getTotalNum() == 1) ? true : false;
    peakfinder_m = std::unique_ptr<PeakFinder> (new PeakFinder(getOutputFN(), rstart_m, rend_m, step_m, singlemode));
}

void Probe::doGoOffline() {
    *gmsg << "* Probe goes offline " << getName() << endl;
    if (online_m && peakfinder_m)
        peakfinder_m->save();
    peakfinder_m.reset(nullptr);
}

void Probe::setStep(double step) {
    step_m = step;
}

double Probe::getStep() const {
    return step_m;
}

bool Probe::doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {
    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);
    // interested in absolute minimum and maximum
    double xmin = std::min(std::abs(rmin(0)), std::abs(rmax(0)));
    double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
    double ymin = std::min(std::abs(rmin(1)), std::abs(rmax(1)));
    double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
    double rbunch_min = std::hypot(xmin, ymin);
    double rbunch_max = std::hypot(xmax, ymax);

    if( rbunch_max > rstart_m - 10.0 && rbunch_min < rend_m + 10.0 ) {
        size_t tempnum = bunch->getLocalNum();

        Vector_t meanP(0.0, 0.0, 0.0);
        for(unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
            for(int d = 0; d < 3; ++d) {
                meanP(d) += bunch->P[i](d);
            }
        }
        reduce(meanP, meanP, OpAddAssign());
        meanP = meanP / Vector_t(bunch->getTotalNum());

        // change probe width depending on step size and angle of particle
        double sk1, sk2, stangle = 0.0;
        if ( B_m == 0.0 ){
            sk1 = meanP(1)/meanP(0);
            if (sk1 == 0.0)
                stangle = 1.0e12;
            else
                stangle = std::abs(1/sk1);
        } else if (meanP(0) == 0.0 ) {
            sk2 = - A_m/B_m;
            if ( sk2 == 0.0 )
                stangle = 1.0e12;
            else
                stangle = std::abs(1/sk2);
        } else {
            sk1 = meanP(1)/meanP(0);
            sk2 = - A_m/B_m;
            stangle = std::abs(( sk1-sk2 )/(1 + sk1*sk2));
        }
        double lstep = (sqrt(1.0-1.0/(1.0+dot(meanP, meanP))) * Physics::c) * tstep*1.0e-6; // [mm]
        double Swidth = lstep / sqrt( 1 + 1/stangle/stangle );
        setGeom(Swidth);
        Vector_t probepoint;

        for(unsigned int i = 0; i < tempnum; ++i) {
            int pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
            if(pflag == 0) continue;
            // calculate closest point at probe -> better to use momentum direction
            // probe: y = -A/B * x - C/B
            // perpendicular line through R: y = B/A * x + R(1) - B/A * R(0)
            // probepoint(0) = (B_m*B_m*bunch->R[i](0) - A_m*B_m*bunch->R[i](1) - A_m*C_m) / (R_m*R_m);
            // probepoint(1) = (A_m*A_m*bunch->R[i](1) - A_m*B_m*bunch->R[i](0) - B_m*C_m) / (R_m*R_m);
            // probepoint(2) = bunch->R[i](2);
            // calculate time correction for probepoint
            // dist1 > 0, right hand, dt > 0; dist1 < 0, left hand, dt < 0
            double dist1 = (A_m*bunch->R[i](0)+B_m*bunch->R[i](1)+C_m)/R_m/1000.0; // [m]
            double k1, k2, tangle = 0.0;
            if ( B_m == 0.0 ){
                k1 = bunch->P[i](1)/bunch->P[i](0);
                if (k1 == 0.0)
                    tangle = 1.0e12;
                else
                    tangle = std::abs(1/k1);
            } else if ( bunch->P[i](0) == 0.0 ) {
                k2 = -A_m/B_m;
                if (k2 == 0.0)
                    tangle = 1.0e12;
                else
                    tangle = std::abs(1/k2);
            } else {
                k1 = bunch->P[i](1)/bunch->P[i](0);
                k2 = -A_m/B_m;
                tangle = std::abs(( k1-k2 )/(1 + k1*k2));
            }
            double dist2 = dist1 * sqrt( 1+1/tangle/tangle );
            double dt = dist2/(sqrt(1.0-1.0/(1.0 + dot(bunch->P[i], bunch->P[i]))) * Physics::c)*1.0e9;

            probepoint = bunch->R[i] + dist2 * 1000.0 * bunch->P[i] / euclidean_norm(bunch->P[i]);

            // peak finder uses millimetre not metre
            peakfinder_m->addParticle(probepoint);

            /*FIXME Issue #45: mm --> m (when OPAL-Cycl uses metre insteas of millimetre,
             * this can be removed.
             */
            probepoint *= 0.001;

            lossDs_m->addParticle(probepoint, bunch->P[i], bunch->ID[i], t+dt, turnnumber);
        }

        peakfinder_m->evaluate(turnnumber);
    }

    // we do not lose particles in the probe
    return false;
}

ElementBase::ElementType Probe::getType() const {
    return PROBE;
}