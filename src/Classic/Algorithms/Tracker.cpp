// ------------------------------------------------------------------------
// $RCSfile: Tracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Tracker
//   The visitor class for tracking a bunch of particles through a beamline
//   using a thin-lens approximation for all elements.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/Tracker.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/Patch.h"
#include "Fields/BMultipoleField.h"

//FIXME Remove headers and dynamic_cast in readOneBunchFromFile
#include "Algorithms/PartBunch.h"
#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

#include <cfloat>
#include <cmath>
#include <limits>

typedef FTps<double, 2> Series2;
typedef FTps<double, 6> Series;

// Class Tracker
// ------------------------------------------------------------------------


Tracker::Tracker(const Beamline &beamline, const PartData &reference,
                 bool backBeam, bool backTrack):
    AbstractTracker(beamline, reference, backBeam, backTrack),
    itsBeamline_m(beamline),
    itsBunch_m(nullptr)
{ }


Tracker::Tracker(const Beamline &beamline,
                 PartBunchBase<double, 3> *bunch,
                 const PartData &reference,
                 bool backBeam, bool backTrack):
    AbstractTracker(beamline, reference, backBeam, backTrack),
    itsBeamline_m(beamline),
    itsBunch_m(bunch)
{}


Tracker::~Tracker()
{}


const PartBunchBase<double, 3> *Tracker::getBunch() const {
    return itsBunch_m;
}


void Tracker::addToBunch(const OpalParticle &part) {
    itsBunch_m->push_back(part);
}


//~ void Tracker::setBunch(const PartBunch &bunch) {
    //~ itsBunch_m = &bunch;
//~ }


void Tracker::visitComponent(const Component &comp) {
    comp.trackBunch(itsBunch_m, itsReference, back_beam, back_track);
}


void Tracker::visitPatch(const Patch &patch) {
    Euclid3D transform = patch.getPatch();
    if(back_path) transform = Inverse(transform);
    applyTransform(transform);
}


void Tracker::visitAlignWrapper(const AlignWrapper &wrap) {
    if(wrap.offset().isIdentity()) {
        wrap.getElement()->accept(*this);
    } else {
        Euclid3D e1 = wrap.getEntranceTransform();
        Euclid3D e2 = wrap.getExitTransform();

        if(back_path) {
            // Tracking right to left.
            applyTransform(Inverse(e2));
            wrap.getElement()->accept(*this);
            applyTransform(Inverse(e1));
        } else {
            // Tracking left to right.
            applyTransform(e1);
            wrap.getElement()->accept(*this);
            applyTransform(e2);
        }
    }
}


void Tracker::applyDrift(double length) {
    double kin = itsReference.getM() / itsReference.getP();
    double refTime = length / itsReference.getBeta();

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);
        if(part.x() != std::numeric_limits<double>::max()) {
            double px = part.px();
            double py = part.py();
            double pt = part.pt() + 1.0;
            double lByPz = length / std::sqrt(pt * pt - px * px - py * py);
            part.x() += px * lByPz;
            part.y() += py * lByPz;
            part.t() += pt * (refTime / std::sqrt(pt * pt + kin * kin) - lByPz);
        }
        itsBunch_m->set_part(part, i);
    }
}


void Tracker::applyThinMultipole
(const BMultipoleField &field, double scale) {
    int order = field.order();

    if(order > 0) {
        for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
            OpalParticle part = itsBunch_m->get_part(i);
            if(part.x() != std::numeric_limits<double>::max()) {
                double x = part.x();
                double y = part.y();
                double kx = + field.normal(order);
                double ky = - field.skew(order);

                int ord = order;
                while(--ord > 0) {
                    double kxt = x * kx - y * ky;
                    double kyt = x * ky + y * kx;
                    kx = kxt + field.normal(ord);
                    ky = kyt - field.skew(ord);
                }
                part.px() -= kx * scale;
                part.py() += ky * scale;
            }
            itsBunch_m->set_part(part, i);
        }
    }
}


void Tracker::applyThinSBend
(const BMultipoleField &field, double scale, double h) {
    Series2 As = buildSBendVectorPotential2D(field, h) * scale;
    Series2 Fx = As.derivative(0);
    Series2 Fy = As.derivative(1);

    // These substitutions work because As depends on x and y only,
    // and not on px or py.
    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);
        FVector<double, 2> z;
        z[0] = part.x();
        z[1] = part.y();
        part.px() -= Fx.evaluate(z);
        part.py() -= Fy.evaluate(z);
        itsBunch_m->set_part(part, i);
    }
}


void Tracker::applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        double kin = itsReference.getM() / itsReference.getP();
        double refTime = refLength / itsReference.getBeta();

        for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
            OpalParticle part = itsBunch_m->get_part(i);
            double px = part.px();
            double py = part.py();
            double pt = part.pt() + 1.0;
            double pz = std::sqrt(pt * pt - px * px - py * py);

            part.px() = euclid.M(0, 0) * px + euclid.M(1, 0) * py + euclid.M(2, 0) * pz;
            part.py() = euclid.M(0, 1) * px + euclid.M(1, 1) * py + euclid.M(2, 1) * pz;
            pz = euclid.M(0, 2) * px + euclid.M(1, 2) * py + euclid.M(2, 2) * pz;

            double x = part.x() - euclid.getX();
            double y = part.y() - euclid.getY();
            double x2 =
                euclid.M(0, 0) * x + euclid.M(1, 0) * y - euclid.M(2, 0) * euclid.getZ();
            double y2 =
                euclid.M(0, 1) * x + euclid.M(1, 1) * y - euclid.M(2, 1) * euclid.getZ();
            double s2 =
                euclid.M(0, 2) * x + euclid.M(1, 2) * y - euclid.M(2, 2) * euclid.getZ();
            double sByPz = s2 / pz;

            double E = std::sqrt(pt * pt + kin * kin);
            part.x() = x2 - sByPz * part.px();
            part.y() = y2 - sByPz * part.py();
            part.t() += pt * (refTime / E + sByPz);
            itsBunch_m->set_part(part, i);
        }
    }
}


Series2 Tracker::
buildMultipoleVectorPotential2D(const BMultipoleField &field) {
    int order = field.order();

    if(order > 0) {
        static const Series2 x = Series2::makeVariable(0);
        static const Series2 y = Series2::makeVariable(1);
        Series2 kx = + field.normal(order) / double(order);
        Series2 ky = - field.skew(order)   / double(order);

        while(order > 1) {
            Series2 kxt = x * kx - y * ky;
            Series2 kyt = x * ky + y * kx;
            order--;
            kx = kxt + field.normal(order) / double(order);
            ky = kyt - field.skew(order)   / double(order);
        }

        Series2 As = x * kx - y * ky;
        As.setTruncOrder(As.getMaxOrder());
        return As;
    } else {
        return Series2(0.0);
    }
}


Series Tracker::
buildMultipoleVectorPotential(const BMultipoleField &field) {
    int order = field.order();

    if(order > 0) {
        static const Series x = Series::makeVariable(X);
        static const Series y = Series::makeVariable(Y);
        Series kx = + field.normal(order) / double(order);
        Series ky = - field.skew(order)   / double(order);

        while(order > 1) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            order--;
            kx = kxt + field.normal(order) / double(order);
            ky = kyt - field.skew(order)   / double(order);
        }

        Series As = x * kx - y * ky;
        As.setTruncOrder(As.getMaxOrder());
        return As;
    } else {
        return Series(0.0);
    }
}


Series2
Tracker::buildSBendVectorPotential2D(const BMultipoleField &field, double h) {
    int order = field.order();
    Series2 As;

    if(order > 0) {
        static const Series2 x = Series2::makeVariable(0);
        static const Series2 y = Series2::makeVariable(1);

        // Construct terms constant and linear in y.
        Series2 Ae = + field.normal(order); // Term even in y.
        Series2 Ao = - field.skew(order);   // Term odd  in y.

        for(int i = order; --i >= 1;) {
            Ae = Ae * x + field.normal(i);
            Ao = Ao * x - field.skew(i);
        }
        Ae.setTruncOrder(Ae.getMaxOrder());
        Ao.setTruncOrder(Ao.getMaxOrder());

        Series2 hx1 = 1. + h * x; // normalized radius
        Ae = + (Ae * hx1).integral(X);
        Ao = - (Ao * hx1);
        // Add terms up to maximum order.
        As = Ae + y * Ao;

        int k = 2;
        if(k <= order) {
            Series2 yp = y * y / 2.0;

            while(true) {
                // Terms even in y.
                Ae = Ae.derivative(0);
                Ae = h * Ae / hx1 - Ae.derivative(0);
                As += Ae * yp;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(0);
                Ao = h * Ao / hx1 - Ao.derivative(0);
                As += Ao * yp;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    return As;
}


Series
Tracker::buildSBendVectorPotential(const BMultipoleField &field, double h) {
    int order = field.order();
    Series As;

    if(order > 0) {
        static const Series x = Series::makeVariable(X);
        static const Series y = Series::makeVariable(Y);

        // Construct terms constant and linear in y.
        Series Ae = + field.normal(order); // Term even in y.
        Series Ao = - field.skew(order);   // Term odd  in y.

        for(int i = order; --i >= 1;) {
            Ae = Ae * x + field.normal(i);
            Ao = Ao * x - field.skew(i);
        }
        Ae.setTruncOrder(Ae.getMaxOrder());
        Ao.setTruncOrder(Ao.getMaxOrder());

        Series hx1 = 1. + h * x; // normalized radius
        Ae = + (Ae * hx1).integral(X);
        Ao = - (Ao * hx1);
        // Add terms up to maximum order.
        As = Ae + y * Ao;

        int k = 2;
        if(k <= order) {
            Series yp = y * y / 2.0;

            while(true) {
                // Terms even in y.
                Ae = Ae.derivative(X);
                Ae = h * Ae / hx1 - Ae.derivative(X);
                As += Ae * yp;
                if(++k > order) break;
                yp *= y / double(k);

                // Terms odd in y.
                Ao = Ao.derivative(X);
                Ao = h * Ao / hx1 - Ao.derivative(X);
                As += Ao * yp;
                if(++k > order) break;
                yp *= y / double(k);
            }
        }
    }

    return As;
}
