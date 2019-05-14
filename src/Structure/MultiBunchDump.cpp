#include "MultiBunchDump.h"

#include <boost/filesystem.hpp>

#include <iomanip>
#include <sstream>

#include "AbstractObjects/OpalData.h"
#include "Utilities/Timer.h"

#include "OPALconfig.h"
#include "Utilities/Util.h"

#include "Ippl.h"

extern Inform *gmsg;

MultiBunchDump::MultiBunchDump()
    : fbase_m(OpalData::getInstance()->getInputBasename())
    , fext_m(".smb")
{
    // in case of non-restart mode we first delete all *.smb files
    this->remove_m();
}


void MultiBunchDump::writeHeader(const std::string& fname) const {
    
    if ( boost::filesystem::exists(fname) ) {
        return;
    }
    
    std::ofstream out;
    this->open_m(out, fname);
    
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");
    
    out << "SDDS1" << std::endl;
    out << "&description\n"
        << indent << "text=\"Multi Bunch Statistics data '" << OpalData::getInstance()->getInputFn()
        << "' " << dateStr << " " << timeStr << "\",\n"
        << indent << "contents=\"multi bunch stat parameters\"\n"
        << "&end\n";
    out << "&parameter\n"
        << indent << "name=processors,\n"
        << indent << "type=long,\n"
        << indent << "description=\"Number of Cores used\"\n"
        << "&end\n";
    out << "&parameter\n"
        << indent << "name=revision,\n"
        << indent << "type=string,\n"
        << indent << "description=\"git revision of opal\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=t,\n"
        << indent << "type=double,\n"
        << indent << "units=ns,\n"
        << indent << "description=\"1 Time\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=azimuth,\n"
        << indent << "type=double,\n"
        << indent << "units=deg,\n"
        << indent << "description=\"2 Azimuth in global coordinates\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=numParticles,\n"
        << indent << "type=long,\n"
        << indent << "units=1,\n"
        << indent << "description=\"3 Number of Macro Particles\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=energy,\n"
        << indent << "type=double,\n"
        << indent << "units=MeV,\n"
        << indent << "description=\"4 Mean Bunch Energy\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=dE,\n"
        << indent << "type=double,\n"
        << indent << "units=MeV,\n"
        << indent << "description=\"5 energy spread of the beam\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_x,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"6 RMS Beamsize in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_y,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"7 RMS Beamsize in y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_s,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"8 RMS Beamsize in s\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_px,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"9 RMS Normalized Momenta in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_py,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"10 RMS Normalized Momenta in y\"\n"
               << "&end\n";
    out << "&column\n"
        << indent << "name=rms_ps,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"11 RMS Normalized Momenta in z\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=emit_x,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"12 Normalized Emittance x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=emit_y,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"13 Normalized Emittance y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=emit_s,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"14 Normalized Emittance s\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=mean_x,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"15 Mean Beam Position in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=mean_y,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"16 Mean Beam Position in y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=mean_s,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"17 Mean Beam Position in s\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=halo_x,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"18 Halo in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=halo_y,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"19 Halo in y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=halo_z,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"20 Halo in z\"\n"
        << "&end\n";
    
    out << "&data\n"
        << indent << "mode=ascii,\n"
        << indent << "no_row_counts=1\n"
        << "&end\n";

    out << Ippl::getNodes() << std::endl;
    out << OPAL_PROJECT_NAME << " "
        << OPAL_PROJECT_VERSION << " git rev. #"
        << Util::getGitRevision() << std::endl;

    close_m(out);
}


void MultiBunchDump::writeData(const beaminfo_t& binfo, short bunch) {

    if ( Ippl::myNode() > 0 )
        return;

    std::stringstream ss;
    ss << fbase_m << "-bunch-"
       << std::setw(4) << std::setfill('0') << bunch
       << fext_m;
    std::string fname = ss.str();

    writeHeader(fname);

    std::ofstream out;
    open_m(out, fname);
    
    unsigned int pwi = 10;
    
    out << binfo.time       << std::setw(pwi) << "\t"
        << binfo.azimuth    << std::setw(pwi) << "\t"
        << binfo.nParticles << std::setw(pwi) << "\t"
        << binfo.ekin       << std::setw(pwi) << "\t"
        << binfo.dEkin      << std::setw(pwi) << "\t"
        << binfo.rrms[0]    << std::setw(pwi) << "\t"
        << binfo.rrms[1]    << std::setw(pwi) << "\t"
        << binfo.rrms[2]    << std::setw(pwi) << "\t"
        << binfo.prms[0]    << std::setw(pwi) << "\t"
        << binfo.prms[1]    << std::setw(pwi) << "\t"
        << binfo.prms[2]    << std::setw(pwi) << "\t"
        << binfo.emit[0]    << std::setw(pwi) << "\t"
        << binfo.emit[1]    << std::setw(pwi) << "\t"
        << binfo.emit[2]    << std::setw(pwi) << "\t"
        << binfo.mean[0]    << std::setw(pwi) << "\t"
        << binfo.mean[1]    << std::setw(pwi) << "\t"
        << binfo.mean[2]    << std::setw(pwi) << "\t"
        << binfo.halo[0]    << std::setw(pwi) << "\t"
        << binfo.halo[1]    << std::setw(pwi) << "\t"
        << binfo.halo[2]    << std::endl;
    
    close_m(out);
}


void MultiBunchDump::open_m(std::ofstream& out,
                            const std::string& fname) const
{
    std::ios::openmode mode = std::ios::out;
    if ( boost::filesystem::exists(fname) ) {
        mode = std::ios::app;
    }
    
    out.open(fname.c_str(), mode);
    out.precision(15);
    out.setf(std::ios::scientific, std::ios::floatfield);
}


void MultiBunchDump::close_m(std::ofstream& out) const {
    if ( out.is_open() ) {
        out.close();
    }
}


void MultiBunchDump::remove_m() const {
    if ( Ippl::myNode() != 0 || OpalData::getInstance()->inRestartRun() ) {
        *gmsg << "* Multi-Bunch Dump (restart mode): Appending "
              << "multi-bunch data to existing smb-files." << endl;
        return;
    }

    namespace fs = boost::filesystem;
    fs::path curDir(fs::current_path());
    std::list<fs::path> files;
    // 1. Mai 2019
    // https://stackoverflow.com/questions/11140483/how-to-get-list-of-files-with-a-specific-extension-in-a-given-folder
    for (const auto & entry : fs::directory_iterator(curDir)) {
        if ( fs::is_regular_file(entry) && entry.path().extension() == fext_m) {
            files.push_back(entry.path());
        }
    }

    // delete *.smb files
    while ( !files.empty() ) {
        fs::path path = files.front();
        fs::remove(path);
        files.pop_front();
    }
    *gmsg << "* Multi-Bunch Dump: Deleted existing multi-bunch smb-files." << endl;
}
