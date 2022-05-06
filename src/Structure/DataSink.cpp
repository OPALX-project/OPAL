//
// Class DataSink
//   This class acts as an observer during the calculation. It generates diagnostic
//   output of the accelerated beam such as statistical beam descriptors of particle
//   positions, momenta, beam phase space (emittance) etc. These are written to file
//   at periodic time steps during the calculation.
//
//   This class also writes the full beam phase space to an H5 file at periodic time
//   steps in the calculation (this period is different from that of the statistical
//   numbers).

//   Class also writes processor load balancing data to file to track parallel
//   calculation efficiency.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Structure/DataSink.h"

#include "OPALconfig.h"

#include "AbstractObjects/OpalData.h"
#include "Fields/Fieldmap.h"
#include "Physics/Units.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/LBalWriter.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "Utility/FieldDebugFunctions.h"

#ifdef __linux__
    #include "MemoryProfiler.h"
#else
    #include "MemoryWriter.h"
#endif

#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

#ifdef ENABLE_AMR
    #include "Structure/GridLBalWriter.h"
#endif

#include <sstream>

DataSink::DataSink()
    : isMultiBunch_m(false)
{
    this->init();
}


DataSink::DataSink(H5PartWrapper *h5wrapper, bool restart, short numBunch)
    : isMultiBunch_m(numBunch > 1)
{
    if (restart && !Options::enableHDF5) {
        throw OpalException("DataSink::DataSink()",
                            "Can not restart when HDF5 is disabled");
    }

    this->init(restart, h5wrapper, numBunch);

    if ( restart )
        rewindLines();
}


DataSink::DataSink(H5PartWrapper *h5wrapper, short numBunch)
    : DataSink(h5wrapper, false, numBunch)
{ }


void DataSink::dumpH5(PartBunchBase<double, 3> *beam, Vector_t FDext[]) const {
    if (!Options::enableHDF5) return;

    h5Writer_m->writePhaseSpace(beam, FDext);
}


int DataSink::dumpH5(PartBunchBase<double, 3> *beam, Vector_t FDext[], double meanEnergy,
                     double refPr, double refPt, double refPz,
                     double refR, double refTheta, double refZ,
                     double azimuth, double elevation, bool local) const
{
    if (!Options::enableHDF5) return -1;

    return h5Writer_m->writePhaseSpace(beam, FDext, meanEnergy, refPr, refPt, refPz,
                                       refR, refTheta, refZ, azimuth, elevation, local);
}


void DataSink::dumpSDDS(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                        const double& azimuth) const
{
    this->dumpSDDS(beam, FDext, losses_t(), azimuth);
}


void DataSink::dumpSDDS(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                        const losses_t &losses, const double& azimuth) const
{
    beam->calcBeamParameters();

    size_t npOutside = 0;
    if (Options::beamHaloBoundary > 0)
        npOutside = beam->calcNumPartsOutside(Options::beamHaloBoundary*beam->get_rrms());

    IpplTimings::startTimer(StatMarkerTimer_m);

    statWriter_m->write(beam, FDext, losses, azimuth, npOutside);

    beam->gatherLoadBalanceStatistics();

    for (size_t i = 0; i < sddsWriter_m.size(); ++i)
        sddsWriter_m[i]->write(beam);

    IpplTimings::stopTimer(StatMarkerTimer_m);
}


void DataSink::storeCavityInformation() {
    if (!Options::enableHDF5) return;

    h5Writer_m->storeCavityInformation();
}


void DataSink::changeH5Wrapper(H5PartWrapper *h5wrapper) {
    if (!Options::enableHDF5) return;

    h5Writer_m->changeH5Wrapper(h5wrapper);
}


void DataSink::writeGeomToVtk(BoundaryGeometry &bg, std::string fn) {
    if (Ippl::myNode() == 0 && Options::enableVTK) {
        bg.writeGeomToVtk (fn);
    }
}


void DataSink::writeImpactStatistics(const PartBunchBase<double, 3> *beam, long long &step, size_t &impact, double &sey_num,
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

        std::unique_ptr<Inform> ofp(new Inform(nullptr, ffn.c_str(), Inform::APPEND, 0));
        Inform &fid = *ofp;
        fid.precision(6);
        fid << std::setiosflags(std::ios::scientific);
        double t = beam->getT() * Units::s2ns;
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


void DataSink::writeMultiBunchStatistics(PartBunchBase<double, 3> *beam,
                                         MultiBunchHandler* mbhandler_p) {
    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    for (short b = 0; b < mbhandler_p->getNumBunch(); ++b) {
        bool isOk = mbhandler_p->calcBunchBeamParameters(beam, b);
        const MultiBunchHandler::beaminfo_t& binfo = mbhandler_p->getBunchInfo(b);
        if (isOk) {
            mbWriter_m[b]->write(beam, binfo);
        }
    }

    for (size_t i = 0; i < sddsWriter_m.size(); ++i)
        sddsWriter_m[i]->write(beam);

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}


void DataSink::setMultiBunchInitialPathLengh(MultiBunchHandler* mbhandler_p) {
    for (short b = 0; b < mbhandler_p->getNumBunch(); ++b) {
        MultiBunchHandler::beaminfo_t& binfo = mbhandler_p->getBunchInfo(b);
        if (mbWriter_m[b]->exists()) {
            binfo.pathlength = mbWriter_m[b]->getLastValue("s");
        } else if (statWriter_m->exists()) {
            binfo.pathlength = statWriter_m->getLastValue("s");
        } else {
            binfo.pathlength = 0.0;
        }
    }
}


void DataSink::rewindLines() {
    unsigned int linesToRewind = 0;

    double spos = h5Writer_m->getLastPosition();
    if (isMultiBunch_m) {
        /* first check if multi-bunch restart
         *
         * first element of vector belongs to first
         * injected bunch in machine --> rewind lines
         * according to that file --> thus rewind in
         * reversed order
         */
        for (std::vector<mbWriter_t>::reverse_iterator rit = mbWriter_m.rbegin();
             rit != mbWriter_m.rend(); ++rit)
        {
            if ((*rit)->exists()) {
                linesToRewind = (*rit)->rewindToSpos(spos);
                (*rit)->replaceVersionString();
            }
        }
    } else if ( statWriter_m->exists() ) {
        // use stat file to get position
        linesToRewind = statWriter_m->rewindToSpos(spos);
        statWriter_m->replaceVersionString();
    }
    h5Writer_m->close();

    // rewind all others
    if ( linesToRewind > 0 ) {
        for (size_t i = 0; i < sddsWriter_m.size(); ++i) {
            sddsWriter_m[i]->rewindLines(linesToRewind);
            sddsWriter_m[i]->replaceVersionString();
        }
    }
}


void DataSink::init(bool restart, H5PartWrapper* h5wrapper, short numBunch) {
    std::string fn = OpalData::getInstance()->getInputBasename();

    lossWrCounter_m = 0;
    StatMarkerTimer_m = IpplTimings::getTimer("Write Stat");

    statWriter_m = statWriter_t(new StatWriter(fn + std::string(".stat"), restart));

    sddsWriter_m.push_back(
        sddsWriter_t(new LBalWriter(fn + std::string(".lbal"), restart))
    );

#ifdef ENABLE_AMR
    if ( Options::amr ) {
        sddsWriter_m.push_back(
            sddsWriter_t(new GridLBalWriter(fn + std::string(".grid"), restart))
        );
    }
#endif

    if ( Options::memoryDump ) {
#ifdef __linux__
        sddsWriter_m.push_back(
            sddsWriter_t(new MemoryProfiler(fn + std::string(".mem"), restart))
        );
#else
        sddsWriter_m.push_back(
            sddsWriter_t(new MemoryWriter(fn + std::string(".mem"), restart))
        );
#endif
    }

    if ( isMultiBunch_m ) {
        initMultiBunchDump(numBunch);
    }

    if ( Options::enableHDF5 ) {
        h5Writer_m = h5Writer_t(new H5Writer(h5wrapper, restart));
    }
}


void DataSink::initMultiBunchDump(short numBunch) {
    bool restart   = OpalData::getInstance()->inRestartRun();
    std::string fn = OpalData::getInstance()->getInputBasename();
    short bunch = mbWriter_m.size();
    while (bunch < numBunch) {
        std::string fname = fn + std::string("-bunch-") +
                            convertToString(bunch, 4) + std::string(".smb");
        mbWriter_m.push_back(
            mbWriter_t(new MultiBunchDump(fname, restart))
        );
        ++bunch;
    }
}
