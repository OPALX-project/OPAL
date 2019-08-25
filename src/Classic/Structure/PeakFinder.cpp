#include "PeakFinder.h"

#include <algorithm>
#include <cmath>
#include <iterator>

#include "AbstractObjects/OpalData.h"
#include "Ippl.h"

PeakFinder::PeakFinder(std::string elem, double min,
                       double max, double binWidth, bool singlemode)
    : element_m(elem)
    , binWidth_m(binWidth)
    , min_m(min)
    , max_m(max)
    , turn_m(0)
    , peakRadius_m(0.0)
    , registered_m(0)
    , singlemode_m(singlemode)
    , first_m(true)
    , finished_m(false)
{
    if (min_m > max_m) {
        std::swap(min_m, max_m);
    }
    // calculate bins, round up so that histogram is large enough (add one for safety)
    nBins_m = static_cast<unsigned int>(std::ceil(( max_m - min_m ) / binWidth_m)) + 1;
}


void PeakFinder::addParticle(const Vector_t& R) {

    double radius = std::hypot(R(0),R(1));
    radius_m.push_back(radius);

    peakRadius_m += radius;
    ++registered_m;
}


void PeakFinder::evaluate(const int& turn) {

    if ( first_m ) {
        turn_m = turn;
        first_m = false;
    }

    if ( turn_m != turn ) {
        finished_m = true;
    }

    bool globFinished = false;

    if ( !singlemode_m )
        allreduce(finished_m, globFinished, 1, std::logical_and<bool>());
    else
        globFinished = finished_m;

    if ( globFinished ) {

        this->computeCentroid_m();

        turn_m = turn;

        // reset
        peakRadius_m = 0.0;
        registered_m = 0;
        finished_m = false;
    }
}


void PeakFinder::save() {

    createHistogram_m();

    // last turn is not yet computed
    this->computeCentroid_m();

    if ( !peaks_m.empty() ) {
        // only rank 0 will go in here

        fn_m   = element_m + std::string(".peaks");
        hist_m = element_m + std::string(".hist");

        INFOMSG("Save " << fn_m << " and " << hist_m << endl);

        if(OpalData::getInstance()->inRestartRun())
            this->append_m();
        else
            this->open_m();

        this->saveASCII_m();

        this->close_m();

    }

    radius_m.clear();
    globHist_m.clear();
}


void PeakFinder::computeCentroid_m() {
    double globPeakRadius = 0.0;
    int globRegister = 0;

    //FIXME inefficient
    if ( !singlemode_m ) {
        reduce(peakRadius_m, globPeakRadius, 1, std::plus<double>());
        reduce(registered_m, globRegister,   1, std::plus<int>());
    } else {
        globPeakRadius = peakRadius_m;
        globRegister = registered_m;
    }

    if ( Ippl::myNode() == 0 ) {
        if ( globRegister > 0 )
            peaks_m.push_back(globPeakRadius / double(globRegister));
    }
}

void PeakFinder::createHistogram_m() {
    /*
     * create local histograms
     */

    globHist_m.resize(nBins_m);
    container_t locHist(nBins_m,0.0);

    double invBinWidth = 1.0 / binWidth_m;
    for(container_t::iterator it = radius_m.begin(); it != radius_m.end(); ++it) {
        int bin = static_cast<int>(std::abs(*it - min_m ) * invBinWidth);
        if (bin < 0 || (unsigned int)bin >= nBins_m) continue; // Probe might save particles outside its boundary
        ++locHist[bin];
    }

    /*
     * create global histograms
     */
    if ( !singlemode_m ) {
        reduce(&(locHist[0]), &(locHist[0]) + locHist.size(),
               &(globHist_m[0]), OpAddAssign());
    } else {
        globHist_m.swap(locHist);
    }

//     reduce(locHist.data(), globHist_m.data(), locHist.size(), std::plus<double>());
}


void PeakFinder::open_m() {
    os_m.open(fn_m.c_str(), std::ios::out);
    hos_m.open(hist_m.c_str(), std::ios::out);
}


void PeakFinder::append_m() {
    os_m.open(fn_m.c_str(), std::ios::app);
    hos_m.open(hist_m.c_str(), std::ios::app);
}


void PeakFinder::close_m() {
    os_m.close();
    hos_m.close();
}


void PeakFinder::saveASCII_m() {
    os_m << "# Peak Radii (mm)" << std::endl;
    for (auto &radius : peaks_m)
        os_m << radius << std::endl;
        
    hos_m << "# Histogram bin counts (min, max, nbins, binsize) "
          << min_m << " mm "
          << max_m << " mm "
          << nBins_m << " "
          << binWidth_m << " mm" << std::endl;
    for (auto binCount : globHist_m)
        hos_m << binCount << std::endl;
}
