//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/DataSink.h"

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Fields/Fieldmap.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPS.h"
#include "Utilities/Timer.h"

#ifdef __linux__
    #include "MemoryProfiler.h"
#endif


#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif



#include "LBalWriter.h"
#include "MemoryWriter.h"

#ifdef ENABLE_AMR
    #include "GridLBalWriter.h"
#endif


#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <queue>
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


void DataSink::dumpH5(EnvelopeBunch &beam, Vector_t FDext[],
                      double sposHead, double sposRef,
                      double sposTail) const
{
    //FIXME https://gitlab.psi.ch/OPAL/src/issues/245
    if (!Options::enableHDF5) return;
    
    h5Writer_m->writePhaseSpace(beam, FDext, sposHead, sposRef, sposTail);
}


void DataSink::dumpSDDS(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                        const double& azimuth) const
{
    this->dumpSDDS(beam, FDext, losses_t(), azimuth);
}


void DataSink::dumpSDDS(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                        const losses_t &losses, const double& azimuth) const
{
    IpplTimings::startTimer(StatMarkerTimer_m);

    statWriter_m->write(beam, FDext, losses, azimuth);

    for (size_t i = 0; i < sddsWriter_m.size(); ++i)
        sddsWriter_m[i]->write(beam);

    IpplTimings::stopTimer(StatMarkerTimer_m);
}


void DataSink::dumpSDDS(EnvelopeBunch &beam, Vector_t FDext[],
                        double sposHead, double sposRef, double sposTail) const
{
    IpplTimings::startTimer(StatMarkerTimer_m);

    statWriter_m->write(beam, FDext, sposHead, sposRef, sposTail);

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
        for(size_t j = 0; j < bg.getNumBFaces(); j++) {
            if(((Geo_mincoords[2] + Geo_hr(2)*i) < bg.TriBarycenters_m[j](2))
               && (bg.TriBarycenters_m[j](2) < (Geo_hr(2)*i + Geo_hr(2) + Geo_mincoords[2]))) {
                prPartLossZ[i] += bg.TriPrPartloss_m[j];
                sePartLossZ[i] += bg.TriSePartloss_m[j];
                fePartLossZ[i] += bg.TriFEPartloss_m[j];
            }

        }
    }
    for(size_t j = 0; j < bg.getNumBFaces(); j++) {
        // fixme: maybe gether particle loss data, i.e., do a reduce() for
        // each triangle in each node befor write to file.
        fidtr << t_step << std::setw(18) << j << std::setw(18)
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


void DataSink::writeGeomToVtk(BoundaryGeometry &bg, std::string fn) {
    if(Ippl::myNode() == 0) {
        bg.writeGeomToVtk (fn);
    }
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
        binfo.pathlength = mbWriter_m[b]->getLastValue("s");
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

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
