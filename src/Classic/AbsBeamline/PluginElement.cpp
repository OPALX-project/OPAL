#include "AbsBeamline/PluginElement.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"

extern Inform *gmsg;

PluginElement::PluginElement():PluginElement("")
{}

PluginElement::PluginElement(const std::string &name):
    Component(name),
    filename_m(""),
    position_m(0.0) {
    setDimensions(0.0, 0.0, 0.0, 0.0);
}

PluginElement::PluginElement(const PluginElement &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m) {
    setDimensions(right.xstart_m, right.xend_m, right.ystart_m, right.yend_m);
}

PluginElement::~PluginElement() {
    if (online_m)
        goOffline();
}

void PluginElement::initialise(PartBunchBase<double, 3> *bunch, double &, double &) {
    initialise(bunch);
}

void PluginElement::initialise(PartBunchBase<double, 3> *bunch) {
    RefPartBunch_m = bunch;
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getOutputFN(), !Options::asciidump));
    // virtual hook
    doInitialise(bunch);
    goOnline(-1e6);
}

void PluginElement::finalise() {
    doFinalise();
    if (online_m)
        goOffline();
}

void PluginElement::goOffline() {
    if (online_m && lossDs_m)
        lossDs_m->save();
    lossDs_m.reset(nullptr);
    doGoOffline();
    online_m = false;
}

bool PluginElement::bends() const {
    return false;
}

bool PluginElement::apply(const size_t &/*i*/, const double &, Vector_t &, Vector_t &) {
    return false;
}

bool PluginElement::applyToReferenceParticle(const Vector_t &, const Vector_t &, const double &, Vector_t &, Vector_t &) {
    return false;
}

void PluginElement::setOutputFN(std::string fn) {
    filename_m = fn;
}

std::string PluginElement::getOutputFN() const {
    if (filename_m == std::string(""))
        return getName();
    else
        return filename_m.substr(0, filename_m.rfind("."));
}

void PluginElement::setDimensions(double xstart, double xend, double ystart, double yend) {
    xstart_m = xstart;
    ystart_m = ystart;
    xend_m   = xend;
    yend_m   = yend;
    rstart_m = std::hypot(xstart, ystart);
    rend_m   = std::hypot(xend,   yend);
    // start position is the one with lowest radius
    if (rstart_m > rend_m) {
        std::swap(xstart_m, xend_m);
        std::swap(ystart_m, yend_m);
        std::swap(rstart_m, rend_m);
    }
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = std::sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;

    // element equation: A*X + B*Y + C = 0
    // point closest to origin https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    double x_close = 0.0;
    if (R_m > 0.0)
        x_close = - A_m * C_m / (R_m * R_m);

    if (x_close > std::min(xstart_m, xend_m) && x_close < std::max(xstart_m, xend_m) )
        rmin_m  = std::abs(C_m) / std::hypot(A_m,B_m);
    else
        rmin_m = rstart_m;
}

void PluginElement::setGeom(const double dist) {

   double slope;
    if (xend_m == xstart_m)
      slope = 1.0e12;
    else
      slope = (yend_m - ystart_m) / (xend_m - xstart_m);

    double coeff2 = sqrt(1 + slope * slope);
    double coeff1 = slope / coeff2;
    double halfdist = dist / 2.0;
    geom_m[0].x = xstart_m - halfdist * coeff1;
    geom_m[0].y = ystart_m + halfdist / coeff2;

    geom_m[1].x = xstart_m + halfdist * coeff1;
    geom_m[1].y = ystart_m - halfdist / coeff2;

    geom_m[2].x = xend_m + halfdist * coeff1;
    geom_m[2].y = yend_m - halfdist  / coeff2;

    geom_m[3].x = xend_m - halfdist * coeff1;
    geom_m[3].y = yend_m + halfdist / coeff2;

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;

    doSetGeom();
}

void PluginElement::changeWidth(PartBunchBase<double, 3> *bunch, int i, const double tstep, const double tangle) {

    constexpr double c_mmtns = Physics::c * 1.0e-6; // m/s --> mm/ns
    double lstep   = euclidean_norm(bunch->P[i]) / Util::getGamma(bunch->P[i]) * c_mmtns * tstep; // [mm]
    double sWidth  = lstep / sqrt( 1 + 1/tangle/tangle );
    setGeom(sWidth);
}

double PluginElement::calculateIncidentAngle(double xp, double yp) const {
    double k1, k2, tangle = 0.0;
    if ( B_m == 0.0 && xp == 0.0) {
        // width is 0.0, keep non-zero
        tangle = 0.1;
    } else if ( B_m == 0.0 ){
        k1 = yp / xp;
        if (k1 == 0.0)
            tangle = 1.0e12;
        else
            tangle = std::abs(1 / k1);
    } else if ( xp == 0.0 ) {
        k2 = - A_m/B_m;
        if ( k2 == 0.0 )
            tangle = 1.0e12;
        else
            tangle = std::abs(1 / k2);
    } else {
        k1 = yp / xp;
        k2 = - A_m / B_m;
        tangle = std::abs(( k1-k2 ) / (1 + k1*k2));
    }
    return tangle;
}

double PluginElement::getXStart() const {
    return xstart_m;
}

double PluginElement::getXEnd() const {
    return xend_m;
}

double PluginElement::getYStart() const {
    return ystart_m;
}

double PluginElement::getYEnd() const {
    return yend_m;
}

bool PluginElement::check(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {
    bool flag = false;
    // check if bunch close
    bool bunchClose = preCheck(bunch);

    if (bunchClose == true) {
        flag = doCheck(bunch, turnnumber, t, tstep); // virtual hook
    }
    // finalise, can have reduce
    flag = finaliseCheck(bunch, flag);

    return flag;
}

void PluginElement::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd   = position_m + 0.005;
}

int PluginElement::checkPoint(const double &x, const double &y) const {
    int    cn = 0;
    for(int i = 0; i < 4; i++) {
        if((   (geom_m[i].y <= y) && (geom_m[i+1].y >  y))
           || ((geom_m[i].y >  y) && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if(x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)
}

void PluginElement::save() {
    OpalData::OPENMODE openMode;
    if (numPassages_m > 0) {
        openMode = OpalData::OPENMODE::APPEND;
    } else {
        openMode = OpalData::getInstance()->getOpenMode();
    }
    lossDs_m->save(1, openMode);
    numPassages_m++;
}

