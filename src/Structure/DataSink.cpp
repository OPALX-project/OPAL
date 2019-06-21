//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/DataSink.h"

#include "OPALconfig.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Fields/Fieldmap.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPS.h"
#include "Utilities/Timer.h"
#include "Util/SDDSParser.h"

#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

#include "H5hut.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <queue>
#include <sstream>

extern Inform *gmsg;

DataSink::DataSink() :
    nMaxBunches_m(1),
    H5call_m(0),
    lossWrCounter_m(0),
    doHDF5_m(true),
    h5wrapper_m(NULL)
#ifdef __linux__
    , memprof_m(new MemoryProfiler())
#endif
{ }

DataSink::DataSink(H5PartWrapper *h5wrapper, int restartStep):
    mode_m(std::ios::out),
    nMaxBunches_m(1),
    H5call_m(0),
    lossWrCounter_m(0),
    h5wrapper_m(h5wrapper)
#ifdef __linux__
    , memprof_m(new MemoryProfiler())
#endif
{
    namespace fs = boost::filesystem;

    doHDF5_m = Options::enableHDF5;
    if (!doHDF5_m) {
        throw OpalException("DataSink::DataSink(int)",
                            "Can not restart when HDF5 is disabled");
    }

    H5PartTimer_m = IpplTimings::getTimer("Write H5-File");
    StatMarkerTimer_m = IpplTimings::getTimer("Write Stat");

    std::string fn = OpalData::getInstance()->getInputBasename();

    statFileName_m = fn + std::string(".stat");
//    mapStatFileName_m = fn + std::string("-map.stat"); //<---


    lBalFileName_m = fn + std::string(".lbal");
    memFileName_m  = fn + std::string(".mem");
#ifdef ENABLE_AMR
    gridLBalFileName_m = fn + std::string(".grid");
#endif

    unsigned int linesToRewind = 0;
    if (fs::exists(statFileName_m)) {
        mode_m = std::ios::app;
        INFOMSG("Appending statistical data to existing data file: " << statFileName_m << endl);
        double spos = h5wrapper->getLastPosition();
        linesToRewind = rewindSDDStoSPos(spos);
        replaceVersionString(statFileName_m);
    } else {
        INFOMSG("Creating new file for statistical data: " << statFileName_m << endl);
    }

//    if (fs::exists(mapStatFileName_m)) //<--
//            mode_m = std::ios::app;
//            INFOMSG("Appending statistical data to existing data file: " << mapStatFileName_m << endl);
//            double spos = h5wrapper->getLastPosition();
//            linesToRewind = rewindSDDStoSPos(spos);
//            replaceVersionString(mapStatFileName_m);
//        } else {
//            INFOMSG("Creating new file for statistical data: " << statFileName_m << endl);
//        }

    if (fs::exists(lBalFileName_m)) {
        INFOMSG("Appending load balance data to existing data file: " << lBalFileName_m << endl);
        rewindLinesLBal(linesToRewind);
    } else {
        INFOMSG("Creating new file for load balance data: " << lBalFileName_m << endl);
    }

    if (Options::memoryDump) {
#ifndef __linux__
        if (fs::exists(memFileName_m)) {
            INFOMSG("Appending memory consumption to existing data file: " << memFileName_m << endl);
            if (Ippl::myNode() == 0) {
                rewindLines(memFileName_m, linesToRewind);
            }
        } else {
            INFOMSG("Creating new file for memory consumption data: " << memFileName_m << endl);
        }
#endif
    }

#ifdef ENABLE_AMR
    if (fs::exists(gridLBalFileName_m)) {
        INFOMSG("Appending grid load balancing to existing data file: " << gridLBalFileName_m << endl);
        if (Ippl::myNode() == 0) {
            rewindLines(gridLBalFileName_m, linesToRewind);
        }
    } else {
        INFOMSG("Creating new file for grid load balancing data: " << gridLBalFileName_m << endl);
    }
#endif

    h5wrapper_m->close();
}

DataSink::DataSink(H5PartWrapper *h5wrapper):
    mode_m(std::ios::out),
    nMaxBunches_m(1),
    H5call_m(0),
    lossWrCounter_m(0),
    h5wrapper_m(h5wrapper)
#ifdef __linux__
    , memprof_m(new MemoryProfiler())
#endif
{
    /// Constructor steps:
    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("Write H5-File");
    StatMarkerTimer_m = IpplTimings::getTimer("Write Stat");

    /// Set file write flags to true. These will be set to false after first
    /// write operation.
    firstWriteH5Surface_m = true;
    /// Define file names.
    std::string fn = OpalData::getInstance()->getInputBasename();
    surfaceLossFileName_m = fn + std::string(".SurfaceLoss.h5");
    statFileName_m = fn + std::string(".stat");
    lBalFileName_m = fn + std::string(".lbal");
    memFileName_m  = fn + std::string(".mem");
#ifdef ENABLE_AMR
    gridLBalFileName_m = fn + std::string(".grid");
#endif

    doHDF5_m = Options::enableHDF5;

    h5wrapper_m->writeHeader();
    h5wrapper_m->close();
}

DataSink::~DataSink() {
    h5wrapper_m = NULL;
}

void DataSink::storeCavityInformation() {
    if (!doHDF5_m) return;

    h5wrapper_m->storeCavityInformation();
}

void DataSink::setMaxNumBunches(int nBunches) {
    nMaxBunches_m = nBunches;
}

void DataSink::writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[]) {

    if (!doHDF5_m) return;

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("B-ref_x", FDext[0](0)),
        std::make_pair("B-ref_z", FDext[0](1)),
        std::make_pair("B-ref_y", FDext[0](2)),
        std::make_pair("E-ref_x", FDext[1](0)),
        std::make_pair("E-ref_z", FDext[1](1)),
        std::make_pair("E-ref_y", FDext[1](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);

    return;
}



int DataSink::writePhaseSpace_cycl(PartBunchBase<double, 3> *beam, Vector_t FDext[], double meanEnergy,
                                   double refPr, double refPt, double refPz,
                                   double refR, double refTheta, double refZ,
                                   double azimuth, double elevation, bool local) {

    if (!doHDF5_m) return -1;
    if (beam->getTotalNum() < 3) return -1; // in single particle mode and tune calculation (2 particles) we do not need h5 data

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("REFPR", refPr),
        std::make_pair("REFPT", refPt),
        std::make_pair("REFPZ", refPz),
        std::make_pair("REFR", refR),
        std::make_pair("REFTHETA", refTheta),
        std::make_pair("REFZ", refZ),
        std::make_pair("AZIMUTH", azimuth),
        std::make_pair("ELEVATION", elevation),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x",  FDext[2](0)),
        std::make_pair("B-ref_z",  FDext[2](1)),
        std::make_pair("B-ref_y",  FDext[2](2)),
        std::make_pair("E-ref_x",  FDext[3](0)),
        std::make_pair("E-ref_z",  FDext[3](1)),
        std::make_pair("E-ref_y",  FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);

    ++ H5call_m;
    return H5call_m - 1;
}

void DataSink::writePhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("sposHead", sposHead),
        std::make_pair("sposRef",  sposRef),
        std::make_pair("sposTail", sposTail),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x",  FDext[2](0)),
        std::make_pair("B-ref_z",  FDext[2](1)),
        std::make_pair("B-ref_y",  FDext[2](2)),
        std::make_pair("E-ref_x",  FDext[3](0)),
        std::make_pair("E-ref_z",  FDext[3](1)),
        std::make_pair("E-ref_y",  FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(&beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::stashPhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    static_cast<H5PartWrapperForPS*>(h5wrapper_m)->stashPhaseSpaceEnvelope(beam,
                                                                           FDext,
                                                                           sposHead,
                                                                           sposRef,
                                                                           sposTail);
    H5call_m++;

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::dumpStashedPhaseSpaceEnvelope() {

    if (!doHDF5_m) return;

    static_cast<H5PartWrapperForPS*>(h5wrapper_m)->dumpStashedPhaseSpaceEnvelope();

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::writeStatData(PartBunchBase<double, 3> *beam,
                             Vector_t FDext[],
                             double E,
                             const double& azimuth)
{
    doWriteStatData(beam, FDext, E, losses_t(), azimuth);
}

void DataSink::writeStatData(PartBunchBase<double, 3> *beam,
                             Vector_t FDext[],
                             const losses_t& losses,
                             const double& azimuth)
{
    doWriteStatData(beam, FDext, beam->get_meanKineticEnergy(), losses, azimuth);
}


void DataSink::doWriteStatData(PartBunchBase<double, 3> *beam,
                               Vector_t FDext[],
                               double Ekin,
                               const losses_t &losses,
                               const double& azimuth)
{

    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    /// Calculate beam statistics and gather load balance statistics.
    beam->calcBeamParameters();
    beam->gatherLoadBalanceStatistics();

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {
        amrbeam->gatherLevelStatistics();
    }
#endif

    size_t npOutside = 0;
    if (Options::beamHaloBoundary>0)
        npOutside = beam->calcNumPartsOutside(Options::beamHaloBoundary*beam->get_rrms());
    // *gmsg << "npOutside 1 = " << npOutside << " beamHaloBoundary= " << Options::beamHaloBoundary << " rrms= " << beam->get_rrms() << endl;

    double  pathLength = 0.0;
    if (OpalData::getInstance()->isInOPALCyclMode())
        pathLength = beam->getLPath();
    else
        pathLength = beam->get_sPos();

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    std::ofstream os_statData;
    std::ofstream os_lBalData;

#ifndef __linux__
    std::ofstream os_memData;
#endif

#ifdef ENABLE_AMR
    std::ofstream os_gridLBalData;
#endif
    double Q = beam->getCharge();

    if ( Options::memoryDump ) {
#ifdef __linux__
        memprof_m->write(beam->getT() * 1e9);
#else
        IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
        memory->sample();
#endif
    }

    if (Ippl::myNode() == 0) {

        open_m(os_statData, statFileName_m);
        open_m(os_lBalData, lBalFileName_m);

#ifndef __linux__
        if ( Options::memoryDump )
            open_m(os_memData, memFileName_m);
#endif

#ifdef ENABLE_AMR
        if ( dynamic_cast<AmrPartBunch*>(beam) != nullptr )
            open_m(os_gridLBalData, gridLBalFileName_m);
#endif
        if (mode_m == std::ios::out) {
            mode_m = std::ios::app;

            writeSDDSHeader(os_statData, losses);
            writeLBalHeader(beam, os_lBalData);

#ifndef __linux__
            if ( Options::memoryDump )
                writeMemoryHeader(os_memData);
#endif

#ifdef ENABLE_AMR
            if ( dynamic_cast<AmrPartBunch*>(beam) != nullptr )
                writeGridLBalHeader(beam, os_gridLBalData);
#endif
        }

        os_statData << beam->getT() * 1e9 << std::setw(pwi) << "\t"         // 1
                    << pathLength << std::setw(pwi) << "\t"                 // 2

                    << beam->getTotalNum() << std::setw(pwi) << "\t"        // 3
                    << Q << std::setw(pwi) << "\t"                          // 4

                    << Ekin << std::setw(pwi) << "\t"                       // 5

                    << beam->get_rrms()(0) << std::setw(pwi) << "\t"        // 6
                    << beam->get_rrms()(1) << std::setw(pwi) << "\t"        // 7
                    << beam->get_rrms()(2) << std::setw(pwi) << "\t"        // 8

                    << beam->get_prms()(0) << std::setw(pwi) << "\t"        // 9
                    << beam->get_prms()(1) << std::setw(pwi) << "\t"        // 10
                    << beam->get_prms()(2) << std::setw(pwi) << "\t"        // 11

                    << beam->get_norm_emit()(0) << std::setw(pwi) << "\t"   // 12
                    << beam->get_norm_emit()(1) << std::setw(pwi) << "\t"   // 13
                    << beam->get_norm_emit()(2) << std::setw(pwi) << "\t"   // 14

                    << beam->get_rmean()(0)  << std::setw(pwi) << "\t"      // 15
                    << beam->get_rmean()(1)  << std::setw(pwi) << "\t"      // 16
                    << beam->get_rmean()(2)  << std::setw(pwi) << "\t"      // 17

                    << beam->RefPartR_m(0) << std::setw(pwi) << "\t"        // 18
                    << beam->RefPartR_m(1) << std::setw(pwi) << "\t"        // 19
                    << beam->RefPartR_m(2) << std::setw(pwi) << "\t"        // 20

                    << beam->RefPartP_m(0) << std::setw(pwi) << "\t"        // 21
                    << beam->RefPartP_m(1) << std::setw(pwi) << "\t"        // 22
                    << beam->RefPartP_m(2) << std::setw(pwi) << "\t"        // 23

                    << beam->get_maxExtent()(0) << std::setw(pwi) << "\t"   // 24
                    << beam->get_maxExtent()(1) << std::setw(pwi) << "\t"   // 25
                    << beam->get_maxExtent()(2) << std::setw(pwi) << "\t"   // 26

            // Write out Courant Snyder parameters.
                    << beam->get_rprms()(0) << std::setw(pwi) << "\t"       // 27
                    << beam->get_rprms()(1) << std::setw(pwi) << "\t"       // 28
                    << beam->get_rprms()(2) << std::setw(pwi) << "\t"       // 29

            // Write out dispersion.
                    << beam->get_Dx() << std::setw(pwi) << "\t"             // 30
                    << beam->get_DDx() << std::setw(pwi) << "\t"            // 31
                    << beam->get_Dy() << std::setw(pwi) << "\t"             // 32
                    << beam->get_DDy() << std::setw(pwi) << "\t"            // 33


            // Write head/reference particle/tail field information.
                    << FDext[0](0) << std::setw(pwi) << "\t"                // 34 B-ref x
                    << FDext[0](1) << std::setw(pwi) << "\t"                // 35 B-ref y
                    << FDext[0](2) << std::setw(pwi) << "\t"                // 36 B-ref z

                    << FDext[1](0) << std::setw(pwi) << "\t"                // 37 E-ref x
                    << FDext[1](1) << std::setw(pwi) << "\t"                // 38 E-ref y
                    << FDext[1](2) << std::setw(pwi) << "\t"                // 39 E-ref z

                    << beam->getdE() << std::setw(pwi) << "\t"              // 40 dE energy spread
                    << beam->getdT() * 1e9 << std::setw(pwi) << "\t"        // 41 dt time step size
                    << npOutside << std::setw(pwi) << "\t";                 // 42 number of particles outside n*sigma

        if(Ippl::getNodes() == 1 && beam->getLocalNum() > 0) {
            os_statData << beam->R[0](0) << std::setw(pwi) << "\t";         // 43 R0_x
            os_statData << beam->R[0](1) << std::setw(pwi) << "\t";         // 44 R0_y
            os_statData << beam->R[0](2) << std::setw(pwi) << "\t";         // 45 R0_z
            os_statData << beam->P[0](0) << std::setw(pwi) << "\t";         // 46 P0_x
            os_statData << beam->P[0](1) << std::setw(pwi) << "\t";         // 47 P0_y
            os_statData << beam->P[0](2) << std::setw(pwi) << "\t";         // 48 P0_z
        }

        if (OpalData::getInstance()->isInOPALCyclMode()) {

            Vector_t halo = beam->get_halo();
            for (int i = 0; i < 3; ++i)
                os_statData << halo(i) << std::setw(pwi) << "\t";

            os_statData << azimuth << std::setw(pwi) << "\t";
        }

        for(size_t i = 0; i < losses.size(); ++ i) {
            os_statData << losses[i].second << std::setw(pwi) << "\t";
        }
        os_statData << std::endl;

        os_statData.close();

        writeLBalData(beam, os_lBalData, pwi);

        os_lBalData.close();

#ifndef __linux__
        if ( Options::memoryDump ) {
            writeMemoryData(beam, os_memData, pwi);
            os_memData.close();
        }
#endif

#ifdef ENABLE_AMR
        if ( dynamic_cast<AmrPartBunch*>(beam) != nullptr ) {
            writeGridLBalData(beam, os_gridLBalData, pwi);
            os_gridLBalData.close();
        }
#endif
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeStatData(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {
    /// Function steps:
    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    /// Calculate beam statistics and gather load balance statistics.
    beam.calcBeamParameters();
    beam.gatherLoadBalanceStatistics();

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    std::ofstream os_statData;
    std::ofstream os_lBalData;

#ifndef __linux__
    std::ofstream os_memData;
#endif

#ifdef ENABLE_AMR
    std::ofstream os_gridLBalData;
#endif
    double en = beam.get_meanKineticEnergy() * 1e-6;

    if ( Options::memoryDump ) {
#ifdef __linux__
        memprof_m->write(beam.getT() * 1e9);
#else
        IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
        memory->sample();
#endif
    }

    if (Ippl::myNode() == 0) {

        open_m(os_statData, statFileName_m);
        open_m(os_lBalData, lBalFileName_m);

#ifndef __linux__
        if ( Options::memoryDump )
            open_m(os_memData, memFileName_m);
#endif

#ifdef ENABLE_AMR
        if ( dynamic_cast<AmrPartBunch*>(&beam) != nullptr )
            open_m(os_gridLBalData, gridLBalFileName_m);
#endif

        if ( mode_m == std::ios::out ) {
            mode_m = std::ios::app;
            writeSDDSHeader(os_statData);
            writeLBalHeader(&beam, os_lBalData);

#ifndef __linux__
            if ( Options::memoryDump )
                writeMemoryHeader(os_memData);
#endif

#ifdef ENABLE_AMR
            if ( dynamic_cast<AmrPartBunch*>(&beam) != nullptr )
                writeGridLBalHeader(&beam, os_gridLBalData);
#endif
        }

        os_statData << beam.getT() << std::setw(pwi) << "\t"                                       // 1
                    << sposRef << std::setw(pwi) << "\t"                                           // 2

                    << beam.getTotalNum() << std::setw(pwi) << "\t"                                // 3
                    << beam.getTotalNum() * beam.getChargePerParticle() << std::setw(pwi) << "\t"  // 4
                    << en << std::setw(pwi) << "\t"                                                // 5

                    << beam.get_rrms()(0) << std::setw(pwi) << "\t"                                // 6
                    << beam.get_rrms()(1) << std::setw(pwi) << "\t"                                // 7
                    << beam.get_rrms()(2) << std::setw(pwi) << "\t"                                // 8

                    << beam.get_prms()(0) << std::setw(pwi) << "\t"                                // 9
                    << beam.get_prms()(1) << std::setw(pwi) << "\t"                                // 10
                    << beam.get_prms()(2) << std::setw(pwi) << "\t"                                // 11

                    << beam.get_norm_emit()(0) << std::setw(pwi) << "\t"                           // 12
                    << beam.get_norm_emit()(1) << std::setw(pwi) << "\t"                           // 13
                    << beam.get_norm_emit()(2) << std::setw(pwi) << "\t"                           // 14

                    << beam.get_rmean()(0)  << std::setw(pwi) << "\t"                              // 15
                    << beam.get_rmean()(1)  << std::setw(pwi) << "\t"                              // 16
                    << beam.get_rmean()(2)  << std::setw(pwi) << "\t"                              // 17

                    << beam.get_maxExtent()(0) << std::setw(pwi) << "\t"                           // 18
                    << beam.get_maxExtent()(1) << std::setw(pwi) << "\t"                           // 19
                    << beam.get_maxExtent()(2) << std::setw(pwi) << "\t"                           // 20

            // Write out Courant Snyder parameters.
                    << 0.0  << std::setw(pwi) << "\t"                                              // 21
                    << 0.0  << std::setw(pwi) << "\t"                                              // 22

                    << 0.0 << std::setw(pwi) << "\t"                                               // 23
                    << 0.0 << std::setw(pwi) << "\t"                                               // 24

            // Write out dispersion.
                    << beam.get_Dx() << std::setw(pwi) << "\t"                                     // 25
                    << beam.get_DDx() << std::setw(pwi) << "\t"                                    // 26
                    << beam.get_Dy() << std::setw(pwi) << "\t"                                     // 27
                    << beam.get_DDy() << std::setw(pwi) << "\t"                                    // 28


            // Write head/reference particle/tail field information.
                    << FDext[2](0) << std::setw(pwi) << "\t"                                       // 29 B-ref x
                    << FDext[2](1) << std::setw(pwi) << "\t"                                       // 30 B-ref y
                    << FDext[2](2) << std::setw(pwi) << "\t"                                       // 31 B-ref z

                    << FDext[3](0) << std::setw(pwi) << "\t"                                       // 32 E-ref x
                    << FDext[3](1) << std::setw(pwi) << "\t"                                       // 33 E-ref y
                    << FDext[3](2) << std::setw(pwi) << "\t"                                       // 34 E-ref z

                    << beam.get_dEdt() << std::setw(pwi) << "\t"                                   // 35 dE energy spread

                    << std::endl;

        os_statData.close();

        writeLBalData(&beam, os_lBalData, pwi);
        os_lBalData.close();

#ifndef __linux__
        if ( Options::memoryDump ) {
            writeMemoryData(&beam, os_memData, pwi);
            os_memData.close();
        }
#endif

#ifdef ENABLE_AMR
        if ( dynamic_cast<AmrPartBunch*>(&beam) != nullptr ) {
            writeGridLBalData(&beam, os_gridLBalData, pwi);
            os_gridLBalData.close();
        }
#endif
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeSDDSHeader(std::ofstream &outputFile) {
    writeSDDSHeader(outputFile,
                    losses_t());
}
void DataSink::writeSDDSHeader(std::ofstream &outputFile,
                               const losses_t &losses) {
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");

    outputFile << "SDDS1" << std::endl;
    outputFile << "&description\n"
               << indent << "text=\"Statistics data '" << OpalData::getInstance()->getInputFn()
               << "' " << dateStr << " " << timeStr << "\",\n"
               << indent << "contents=\"stat parameters\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=processors,\n"
               << indent << "type=long,\n"
               << indent << "description=\"Number of Cores used\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=revision,\n"
               << indent << "type=string,\n"
               << indent << "description=\"git revision of opal\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=flavor,\n"
               << indent << "type=string,\n"
               << indent << "description=\"OPAL flavor that wrote file\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=t,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"1 Time\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=s,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"2 Path length\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=numParticles,\n"
               << indent << "type=long,\n"
               << indent << "units=1,\n"
               << indent << "description=\"3 Number of Macro Particles\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=charge,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"4 Bunch Charge\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=energy,\n"
               << indent << "type=double,\n"
               << indent << "units=MeV,\n"
               << indent << "description=\"5 Mean Bunch Energy\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=rms_x,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"6 RMS Beamsize in x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=rms_y,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"7 RMS Beamsize in y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=rms_s,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"8 RMS Beamsize in s\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=rms_px,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"9 RMS Normalized Momenta in x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=rms_py,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"10 RMS Normalized Momenta in y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=rms_ps,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"11 RMS Normalized Momenta in s\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=emit_x,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"12 Normalized Emittance x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=emit_y,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"13 Normalized Emittance y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=emit_s,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"14 Normalized Emittance s\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=mean_x,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"15 Mean Beam Position in x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=mean_y,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"16 Mean Beam Position in y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=mean_s,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"17 Mean Beam Position in s\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ref_x,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"18 x coordinate of reference particle in lab cs\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ref_y,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"19 y coordinate of reference particle in lab cs\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ref_z,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"20 z coordinate of reference particle in lab cs\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ref_px,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"21 x momentum of reference particle in lab cs\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ref_py,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"22 y momentum of reference particle in lab cs\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ref_pz,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"23 z momentum of reference particle in lab cs\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=max_x,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"24 Max Beamsize in x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=max_y,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"25 Max Beamsize in y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=max_s,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"26 Max Beamsize in s\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=xpx,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"27 Correlation xpx\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=ypy,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"28 Correlation ypy\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=zpz,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"29 Correlation zpz\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Dx,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"30 Dispersion in x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=DDx,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"31 Derivative of dispersion in x\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Dy,\n"
               << indent << "type=double,\n"
               << indent << "units=m,\n"
               << indent << "description=\"32 Dispersion in y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=DDy,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"33 Derivative of dispersion in y\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Bx_ref,\n"
               << indent << "type=double,\n"
               << indent << "units=T,\n"
               << indent << "description=\"34 Bx-Field component of ref particle\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=By_ref,\n"
               << indent << "type=double,\n"
               << indent << "units=T,\n"
               << indent << "description=\"35 By-Field component of ref particle\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Bz_ref,\n"
               << indent << "type=double,\n"
               << indent << "units=T,\n"
               << indent << "description=\"36 Bz-Field component of ref particle\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Ex_ref,\n"
               << indent << "type=double,\n"
               << indent << "units=MV/m,\n"
               << indent << "description=\"37 Ex-Field component of ref particle\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Ey_ref,\n"
               << indent << "type=double,\n"
               << indent << "units=MV/m,\n"
               << indent << "description=\"38 Ey-Field component of ref particle\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=Ez_ref,\n"
               << indent << "type=double,\n"
               << indent << "units=MV/m,\n"
               << indent << "description=\"39 Ez-Field component of ref particle\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=dE,\n"
               << indent << "type=double,\n"
               << indent << "units=MeV,\n"
               << indent << "description=\"40 energy spread of the beam\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=dt,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"41 time step size\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=partsOutside,\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"42 outside n*sigma of the beam\"\n"
               << "&end\n";

    unsigned int columnStart = 43;
    if(Ippl::getNodes() == 1) {
        outputFile << "&column\n"
                   << indent << "name=R0_x,\n"
                   << indent << "type=double,\n"
                   << indent << "units=m,\n"
                   << indent << "description=\"43 R0 Particle position in x\"\n"
                   << "&end\n";
        outputFile << "&column\n"
                   << indent << "name=R0_y,\n"
                   << indent << "type=double,\n"
                   << indent << "units=m,\n"
                   << indent << "description=\"44 R0 Particle position in y\"\n"
                   << "&end\n";
        outputFile << "&column\n"
                   << indent << "name=R0_s,\n"
                   << indent << "type=double,\n"
                   << indent << "units=m,\n"
                   << indent << "description=\"45 R0 Particle position in z\"\n"
                   << "&end\n";
        outputFile << "&column\n"
                   << indent << "name=P0_x,\n"
                   << indent << "type=double,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"46 R0 Particle momentum in x\"\n"
                   << "&end\n";
        outputFile << "&column\n"
                   << indent << "name=P0_y,\n"
                   << indent << "type=double,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"47 R0 Particle momentum in y\"\n"
                   << "&end\n";
        outputFile << "&column\n"
                   << indent << "name=P0_s,\n"
                   << indent << "type=double,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"48 R0 Particle momentum in z\"\n"
                   << "&end\n";
        columnStart = 49;
    }

    if (OpalData::getInstance()->isInOPALCyclMode()) {
        char dir[] = { 'x', 'y', 'z' };

        for (int i = 0; i < 3; ++i) {
            std::stringstream ss;
            ss << "&column\n" << indent << "name=halo_" << dir[i] << ",\n"
               << indent << "type=double,\n"
               << indent << "units=1,\n"
               << indent << "description=\"" << columnStart++ << " Halo in "
               << dir[i] << "\"\n"
               << "&end\n";
               outputFile << ss.str();
        }
        std::stringstream ss;
        ss << "&column\n" << indent << "name=azimuth,\n"
           << indent << "type=double,\n"
           << indent << "units=deg,\n"
           << indent << "description=\"" << columnStart++ << " Azimuth in "
           << "global coordinates\"\n"
           << "&end\n";
        outputFile << ss.str();
    }

    for (size_t i = 0; i < losses.size(); ++ i) {
        outputFile << "&column\n"
                   << indent << "name=" << losses[i].first << ",\n"
                   << indent << "type=long,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"" << columnStart ++ << "Number of lost particles in element\"\n"
                   << "&end\n";
    }

//    std::string method = Util::toUpper(Attributes::getString(itsAttr[METHOD])); //<--
//    if (method == "THICK") {
//        outputFile << "&column\n"
//               << indent << "name=mapD_x,\n"
//               << indent << "type=double,\n"
//               << indent << "units=1,\n"
//               << indent << "description=\"49 Dx dispersion of beamline"\n"
//               << "&end\n";
//        outputFile << "&column\n"
//               << indent << "name=mapD_y,\n"
//               << indent << "type=double,\n"
//               << indent << "units=1,\n"
//               << indent << "description=\"50 Dy dispersion of beamline"\n"
//               << "&end\n";
//        columnStart = 51;
//    }

    outputFile << "&data\n"
               << indent << "mode=ascii,\n"
               << indent << "no_row_counts=1\n"
               << "&end\n";

    outputFile << Ippl::getNodes() << std::endl;
    outputFile << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl;
    outputFile << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                   (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
}


void DataSink::writePartlossZASCII(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg, std::string fn) {

    size_t temp = lossWrCounter_m ;

    std::string ffn = fn + convertToString(temp) + std::string("Z.dat");
    std::unique_ptr<Inform> ofp(new Inform(NULL, ffn.c_str(), Inform::OVERWRITE, 0));
    Inform &fid = *ofp;
    setInform(fid);
    fid.precision(6);

    std::string ftrn =  fn + std::string("triangle") + convertToString(temp) + std::string(".dat");
    std::unique_ptr<Inform> oftr(new Inform(NULL, ftrn.c_str(), Inform::OVERWRITE, 0));
    Inform &fidtr = *oftr;
    setInform(fidtr);
    fidtr.precision(6);

    Vector_t Geo_nr = bg.getnr();
    Vector_t Geo_hr = bg.gethr();
    Vector_t Geo_mincoords = bg.getmincoords();
    double t = beam->getT();
    double t_step = t * 1.0e9;
    double* prPartLossZ = new double[bg.getnr() (2)];
    double* sePartLossZ = new double[bg.getnr() (2)];
    double* fePartLossZ = new double[bg.getnr() (2)];
    fidtr << "# Time/ns" << std::setw(18) << "Triangle_ID" << std::setw(18)
          << "Xcoordinates (m)" << std::setw(18)
          << "Ycoordinates (m)" << std::setw(18)
          << "Zcoordinates (m)" << std::setw(18)
          << "Primary part. charge (C)" << std::setw(40)
          << "Field emit. part. charge (C)" << std::setw(40)
          << "Secondary emit. part. charge (C)" << std::setw(40) << endl;
    for(int i = 0; i < Geo_nr(2) ; i++) {
        prPartLossZ[i] = 0;
        sePartLossZ[i] = 0;
        fePartLossZ[i] = 0;
        for(int j = 0; j < bg.getNumBFaces(); j++) {
            if(((Geo_mincoords[2] + Geo_hr(2)*i) < bg.TriBarycenters_m[j](2))
               && (bg.TriBarycenters_m[j](2) < (Geo_hr(2)*i + Geo_hr(2) + Geo_mincoords[2]))) {
                prPartLossZ[i] += bg.TriPrPartloss_m[j];
                sePartLossZ[i] += bg.TriSePartloss_m[j];
                fePartLossZ[i] += bg.TriFEPartloss_m[j];
            }

        }
    }
    for(int j = 0; j < bg.getNumBFaces(); j++) {
        fidtr << t_step << std::setw(18) << j << std::setw(18)// fixme: maybe gether particle loss data, i.e., do a reduce() for each triangle in each node befor write to file.
              << bg.TriBarycenters_m[j](0) << std::setw(18)
              << bg.TriBarycenters_m[j](1) << std::setw(18)
              << bg.TriBarycenters_m[j](2) <<  std::setw(40)
              << -bg.TriPrPartloss_m[j] << std::setw(40)
              << -bg.TriFEPartloss_m[j] <<  std::setw(40)
              << -bg.TriSePartloss_m[j] << endl;
    }
    fid << "# Delta_Z/m" << std::setw(18)
        << "Zcoordinates (m)" << std::setw(18)
        << "Primary part. charge (C)" << std::setw(40)
        << "Field emit. part. charge (C)" << std::setw(40)
        << "Secondary emit. part. charge (C)" << std::setw(40) << "t" << endl;


    for(int i = 0; i < Geo_nr(2) ; i++) {
        double primaryPLoss = -prPartLossZ[i];
        double secondaryPLoss = -sePartLossZ[i];
        double fieldemissionPLoss = -fePartLossZ[i];
        reduce(primaryPLoss, primaryPLoss, OpAddAssign());
        reduce(secondaryPLoss, secondaryPLoss, OpAddAssign());
        reduce(fieldemissionPLoss, fieldemissionPLoss, OpAddAssign());
        fid << Geo_hr(2) << std::setw(18)
            << Geo_mincoords[2] + Geo_hr(2)*i << std::setw(18)
            << primaryPLoss << std::setw(40)
            << fieldemissionPLoss << std::setw(40)
            << secondaryPLoss << std::setw(40) << t << endl;
    }
    lossWrCounter_m++;
    delete[] prPartLossZ;
    delete[] sePartLossZ;
    delete[] fePartLossZ;
}

void DataSink::writeSurfaceInteraction(PartBunchBase<double, 3> *beam, long long &step, BoundaryGeometry &bg, std::string fn) {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);
    if(firstWriteH5Surface_m) {
        firstWriteH5Surface_m = false;

        h5_prop_t props = H5CreateFileProp ();
        MPI_Comm comm = Ippl::getComm();
        H5SetPropFileMPIOCollective (props, &comm);
        H5fileS_m = H5OpenFile (surfaceLossFileName_m.c_str(), H5_O_WRONLY, props);
        if(H5fileS_m == (h5_file_t)H5_ERR) {
            throw OpalException("DataSink::writeSurfaceInteraction",
                                "failed to open h5 file '" + surfaceLossFileName_m + "' for surface loss");
        }
        H5CloseProp (props);

    }
    int nTot = bg.getNumBFaces();

    int N_mean = static_cast<int>(floor(nTot / Ippl::getNodes()));
    int N_extra = static_cast<int>(nTot - N_mean * Ippl::getNodes());
    int pc = 0;
    int count = 0;
    if(Ippl::myNode() == 0) {
        N_mean += N_extra;
    }
    std::unique_ptr<char[]> varray(new char[(N_mean)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());


    rc = H5SetStep(H5fileS_m, step);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5PartSetNumParticles(H5fileS_m, N_mean);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    double    qi = beam->getChargePerParticle();
    rc = H5WriteStepAttribFloat64(H5fileS_m, "qi", &qi, 1);

    std::unique_ptr<double[]> tmploss(new double[nTot]);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;

    // :NOTE: may be removed if we have parallelized the geometry .
    reduce(bg.TriPrPartloss_m, bg.TriPrPartloss_m + nTot, tmploss.get(), OpAddAssign());

    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0) {
                    //farray[count] =  bg.TriPrPartloss_m[Ippl::myNode()*N_mean+count];
                    size_t idx = pc * N_mean + count + N_extra;
                    if (((bg.TriBGphysicstag_m[idx] & (BGphysics::Absorption)) == (BGphysics::Absorption)) &&
                        ((bg.TriBGphysicstag_m[idx] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) &&
                        ((bg.TriBGphysicstag_m[idx] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if (((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) &&
                        ((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) &&
                        ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "PrimaryLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;
    reduce(bg.TriSePartloss_m, bg.TriSePartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we parallelize the geometry as well.
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0) {
                    size_t idx = pc * N_mean + count + N_extra;
                    if (((bg.TriBGphysicstag_m[idx] & (BGphysics::Absorption)) == (BGphysics::Absorption)) &&
                        ((bg.TriBGphysicstag_m[idx] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) &&
                        ((bg.TriBGphysicstag_m[idx] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if (((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) &&
                        ((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) &&
                        ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "SecondaryLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;

    reduce(bg.TriFEPartloss_m, bg.TriFEPartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we parallelize the geometry as well.
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {

                if(pc != 0) {
                    size_t idx = pc * N_mean + count + N_extra;
                    if (((bg.TriBGphysicstag_m[idx] & (BGphysics::Absorption)) == (BGphysics::Absorption)) &&
                        ((bg.TriBGphysicstag_m[idx] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) &&
                        ((bg.TriBGphysicstag_m[idx] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if (((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) &&
                        ((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) &&
                        ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }

            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "FNEmissionLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0)
                    larray[count] =  Ippl::myNode() * N_mean + count + N_extra; // node 0 will be 0*N_mean+count+N_extra also correct.
                else
                    larray[count] = count;
                count ++;
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataInt64(H5fileS_m, "TriangleID", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);

}

void DataSink::writeImpactStatistics(PartBunchBase<double, 3> *beam, long long &step, size_t &impact, double &sey_num,
                                     size_t numberOfFieldEmittedParticles, bool nEmissionMode, std::string fn) {

    double charge = 0.0;
    size_t Npart = 0;
    double Npart_d = 0.0;
    if(!nEmissionMode) {
        charge = -1.0 * beam->getCharge();
        //reduce(charge, charge, OpAddAssign());
        Npart_d = -1.0 * charge / beam->getChargePerParticle();
    } else {
        Npart = beam->getTotalNum();
    }
    if(Ippl::myNode() == 0) {
        std::string ffn = fn + std::string(".dat");

        std::unique_ptr<Inform> ofp(new Inform(NULL, ffn.c_str(), Inform::APPEND, 0));
        Inform &fid = *ofp;
        setInform(fid);

        fid.precision(6);
        fid << std::setiosflags(std::ios::scientific);
        double t = beam->getT() * 1.0e9;
        if(!nEmissionMode) {

            if(step == 0) {
                fid << "#Time/ns"  << std::setw(18) << "#Geometry impacts" << std::setw(18) << "tot_sey" << std::setw(18)
                    << "TotalCharge" << std::setw(18) << "PartNum" << " numberOfFieldEmittedParticles " << endl;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num << std::setw(18) << charge
                << std::setw(18) << Npart_d << std::setw(18) << numberOfFieldEmittedParticles << endl;
        } else {

            if(step == 0) {
                fid << "#Time/ns"  << std::setw(18) << "#Geometry impacts" << std::setw(18) << "tot_sey" << std::setw(18)
                    << "ParticleNumber" << " numberOfFieldEmittedParticles " << endl;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num
                << std::setw(18) << double(Npart) << std::setw(18) << numberOfFieldEmittedParticles << endl;
        }
    }
}

void DataSink::writeGeomToVtk(BoundaryGeometry &bg, std::string fn) {
    if(Ippl::myNode() == 0) {
        bg.writeGeomToVtk (fn);
    }
}

/** \brief Find out which if we write HDF5 or not
 *
 *
 *
 */
bool DataSink::doHDF5() {
    return doHDF5_m;
}


/** \brief
 *  delete the last 'numberOfLines' lines of the file 'fileName'
 */
void DataSink::rewindLines(const std::string &fileName, size_t numberOfLines) const {
    if (numberOfLines == 0) return;

    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;

    fs.open (fileName.c_str(), std::fstream::in);

    if (!fs.is_open()) return;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fileName.c_str(), std::fstream::out);

    if (!fs.is_open()) return;

    while (allLines.size() > numberOfLines) {
        fs << allLines.front() << "\n";
        allLines.pop();
    }
    fs.close();
}

void DataSink::replaceVersionString(const std::string &fileName) const {

    if (Ippl::myNode() == 0) {
        std::string versionFile;
        SDDS::SDDSParser parser(fileName);
        parser.run();
        parser.getParameterValue("revision", versionFile);

        std::string line;
        std::queue<std::string> allLines;
        std::fstream fs;

        fs.open (fileName.c_str(), std::fstream::in);

        if (!fs.is_open()) return;

        while (getline(fs, line)) {
            allLines.push(line);
        }
        fs.close();


        fs.open (fileName.c_str(), std::fstream::out);

        if (!fs.is_open()) return;

        while (allLines.size() > 0) {
            line = allLines.front();

            if (line != versionFile) {
                fs << line << "\n";
            } else {
                fs << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << "\n";
            }

            allLines.pop();
        }

        fs.close();
    }
}


void DataSink::open_m(std::ofstream& os, const std::string& fileName) const {
    os.open(fileName.c_str(), mode_m);
    os.precision(15);
    os.setf(std::ios::scientific, std::ios::floatfield);
}


void DataSink::writeLBalHeader(PartBunchBase<double, 3> *beam,
                               std::ofstream &outputFile)
{
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");

    outputFile << "SDDS1" << std::endl;
    outputFile << "&description\n"
               << indent << "text=\"Processor statistics '"
               << OpalData::getInstance()->getInputFn() << "' "
               << dateStr << "" << timeStr << "\",\n"
               << indent << "contents=\"stat parameters\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=processors,\n"
               << indent << "type=long,\n"
               << indent << "description=\"Number of Cores used\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=revision,\n"
               << indent << "type=string,\n"
               << indent << "description=\"git revision of opal\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=flavor,\n"
               << indent << "type=string,\n"
               << indent << "description=\"OPAL flavor that wrote file\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=t,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"1 Time\"\n"
               << "&end\n";

    unsigned int columnStart = 2;


    for (int p = 0; p < Ippl::getNodes(); ++p) {
        outputFile << "&column\n"
                   << indent << "name=processor-" << p << ",\n"
                   << indent << "type=long,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"" << columnStart
                   << " Number of particles of processor " << p << "\"\n"
                   << "&end\n";
        ++columnStart;
    }

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {

        int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

        for (int lev = 0; lev < nLevel; ++lev) {
            outputFile << "&column\n"
                   << indent << "name=level-" << lev << ",\n"
                   << indent << "type=long,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"" << columnStart
                   << " Number of particles at level " << lev << "\"\n"
                   << "&end\n";
            ++columnStart;
        }
    }
#endif

    outputFile << "&data\n"
               << indent << "mode=ascii,\n"
               << indent << "no_row_counts=1\n"
               << "&end\n";

    outputFile << Ippl::getNodes() << std::endl;
    outputFile << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl;
    outputFile << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                   (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;


}

void DataSink::writeLBalData(PartBunchBase<double, 3> *beam,
                             std::ofstream &os_lBalData,
                             unsigned int pwi)
{
    os_lBalData << beam->getT() * 1e9 << std::setw(pwi) << "\t";     // 1

    size_t nProcs = Ippl::getNodes();
    for (size_t p = 0; p < nProcs; ++ p) {
        os_lBalData << beam->getLoadBalance(p)  << std::setw(pwi);

        if ( p + 1 < nProcs )
            os_lBalData << "\t";

    }

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {
        os_lBalData << "\t";
        int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;
        for (int lev = 0; lev < nLevel; ++lev) {
            os_lBalData << amrbeam->getLevelStatistics(lev) << std::setw(pwi);

            if ( lev < nLevel - 1 )
                os_lBalData << "\t";
        }
    }
#endif
    os_lBalData << std::endl;
}


void DataSink::writeMemoryHeader(std::ofstream &outputFile)
{
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    outputFile << "SDDS1" << std::endl;
    outputFile << "&description\n"
               << indent << "text=\"Memory statistics '"
               << OpalData::getInstance()->getInputFn() << "' "
               << dateStr << "" << timeStr << "\",\n"
               << indent << "contents=\"stat parameters\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=processors,\n"
               << indent << "type=long,\n"
               << indent << "description=\"Number of Cores used\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=revision,\n"
               << indent << "type=string,\n"
               << indent << "description=\"git revision of opal\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=flavor,\n"
               << indent << "type=string,\n"
               << indent << "description=\"OPAL flavor that wrote file\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=t,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"1 Time\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=memory,\n"
               << indent << "type=double,\n"
               << indent << "units=" + memory->getUnit() + ",\n"
               << indent << "description=\"2 Total Memory\"\n"
               << "&end\n";

    unsigned int columnStart = 3;

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        outputFile << "&column\n"
                   << indent << "name=processor-" << p << ",\n"
                   << indent << "type=double,\n"
                   << indent << "units=" + memory->getUnit() + ",\n"
                   << indent << "description=\"" << columnStart
                   << " Memory per processor " << p << "\"\n"
                   << "&end\n";
        ++columnStart;
    }

    outputFile << "&data\n"
               << indent << "mode=ascii,\n"
               << indent << "no_row_counts=1\n"
               << "&end\n";

    outputFile << Ippl::getNodes() << std::endl;
    outputFile << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl;
    outputFile << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                   (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
}


void DataSink::writeMemoryData(PartBunchBase<double, 3> *beam,
                               std::ofstream &os_memData,
                               unsigned int pwi)
{
    os_memData << beam->getT() * 1e9 << std::setw(pwi) << "\t";     // 1

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    int nProcs = Ippl::getNodes();
    double total = 0.0;
    for (int p = 0; p < nProcs; ++p) {
        total += memory->getMemoryUsage(p);
    }

    os_memData << total << std::setw(pwi) << "\t";

    for (int p = 0; p < nProcs; p++) {
        os_memData << memory->getMemoryUsage(p)  << std::setw(pwi);

        if ( p < nProcs - 1 )
            os_memData << "\t";

    }
    os_memData << std::endl;
}

#ifdef ENABLE_AMR
void DataSink::writeGridLBalHeader(PartBunchBase<double, 3> *beam,
                                   std::ofstream &outputFile)
{
    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("DataSink::writeGridLBalHeader()",
                            "Can not write grid load balancing for non-AMR runs.");

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");

    outputFile << "SDDS1" << std::endl;
    outputFile << "&description\n"
               << indent << "text=\"Grid load balancing statistics '"
               << OpalData::getInstance()->getInputFn() << "' "
               << dateStr << "" << timeStr << "\",\n"
               << indent << "contents=\"stat parameters\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=processors,\n"
               << indent << "type=long,\n"
               << indent << "description=\"Number of Cores used\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=revision,\n"
               << indent << "type=string,\n"
               << indent << "description=\"git revision of opal\"\n"
               << "&end\n";
    outputFile << "&parameter\n"
               << indent << "name=flavor,\n"
               << indent << "type=string,\n"
               << indent << "description=\"OPAL flavor that wrote file\"\n"
               << "&end\n";
    outputFile << "&column\n"
               << indent << "name=t,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"1 Time\"\n"
               << "&end\n";

    unsigned int columnStart = 2;

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

    for (int lev = 0; lev < nLevel; ++lev) {
        outputFile << "&column\n"
                   << indent << "name=level-" << lev << ",\n"
                   << indent << "type=long,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"" << columnStart
                   << " Number of boxes at level " << lev << "\"\n"
                   << "&end\n";
        ++columnStart;
    }


    for (int p = 0; p < Ippl::getNodes(); ++p) {
        outputFile << "&column\n"
                   << indent << "name=processor-" << p << ",\n"
                   << indent << "type=long,\n"
                   << indent << "units=1,\n"
                   << indent << "description=\"" << columnStart
                   << " Number of grid points per processor " << p << "\"\n"
                   << "&end\n";
        ++columnStart;
    }

    outputFile << "&data\n"
               << indent << "mode=ascii,\n"
               << indent << "no_row_counts=1\n"
               << "&end\n";

    outputFile << Ippl::getNodes() << std::endl;
    outputFile << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl;
    outputFile << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                   (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
}


void DataSink::writeGridLBalData(PartBunchBase<double, 3> *beam,
                                 std::ofstream &os_gridLBalData,
                                 unsigned int pwi)
{
    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("DataSink::writeGridLBalData()",
                            "Can not write grid load balancing for non-AMR runs.");

    os_gridLBalData << amrbeam->getT() * 1e9 << std::setw(pwi) << "\t";     // 1

    std::map<int, long> gridPtsPerCore;

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;
    std::vector<int> gridsPerLevel;

    amrbeam->getAmrObject()->getGridStatistics(gridPtsPerCore, gridsPerLevel);

    os_gridLBalData << "\t";
    for (int lev = 0; lev < nLevel; ++lev) {
        os_gridLBalData << gridsPerLevel[lev] << std::setw(pwi) << "\t";
    }

    int nProcs = Ippl::getNodes();
    for (int p = 0; p < nProcs; ++p) {
        os_gridLBalData << gridPtsPerCore[p] << std::setw(pwi);

        if ( p < nProcs - 1 )
            os_gridLBalData << "\t";
    }
    os_gridLBalData << std::endl;
}


bool DataSink::writeAmrStatistics(PartBunchBase<double, 3> *beam) {

    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        return false;

    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    beam->gatherLoadBalanceStatistics();

    amrbeam->gatherLevelStatistics();

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    std::ofstream os_lBalData;
#ifndef __linux__
    std::ofstream os_memData;
#endif
    std::ofstream os_gridLBalData;

    if ( Options::memoryDump ) {
#ifdef __linux__
        memprof_m->write(beam->getT() * 1e9);
#else
        IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
        memory->sample();
#endif
    }

    if (Ippl::myNode() == 0) {
        open_m(os_lBalData, lBalFileName_m);

#ifndef __linux__
        if ( Options::memoryDump )
            open_m(os_memData, memFileName_m);
#endif

        open_m(os_gridLBalData, gridLBalFileName_m);

        if (mode_m == std::ios::out) {
            mode_m = std::ios::app;

            writeLBalHeader(beam, os_lBalData);

#ifndef __linux__
            if ( Options::memoryDump )
                writeMemoryHeader(os_memData);
#endif

            writeGridLBalHeader(beam, os_gridLBalData);
        }

        writeLBalData(beam, os_lBalData, pwi);

        os_lBalData.close();

#ifndef __linux__
        if ( Options::memoryDump ) {
            writeMemoryData(beam, os_memData, pwi);
            os_memData.close();
        }
#endif

        writeGridLBalData(beam, os_gridLBalData, pwi);
        os_gridLBalData.close();
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
    
    return true;
}


void DataSink::noAmrDump(PartBunchBase<double, 3> *beam) {
#ifdef __linux__
    if ( Options::memoryDump ) {
        memprof_m->write(beam->getT() * 1e9);
    }
#else
    if ( Ippl::myNode() != 0 ) {
        return;
    }

    std::ofstream os_memData;
    if ( Options::memoryDump ) {
        open_m(os_memData, memFileName_m);
    }

    std::ofstream os_lBalData;
    open_m(os_lBalData, lBalFileName_m);

    if (mode_m == std::ios::out) {
        mode_m = std::ios::app;

        writeLBalHeader(os_lBalData);

        if ( Options::memoryDump )
            writeMemoryHeader(os_memData);
    }

    unsigned int pwi = 10;

    if ( Options::memoryDump ) {
        writeMemoryData(beam, os_memData, pwi);
        os_memData.close();
    }

    writeLBalData(beam, os_lBalData, pwi);
    os_lBalData.close();
#endif
}

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: