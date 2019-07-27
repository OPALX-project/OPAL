#include "MultiBunchDump.h"

#include <sstream>

#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Ippl.h"

extern Inform *gmsg;

MultiBunchDump::MultiBunchDump(const std::string& fname, bool restart)
    : StatBaseWriter(fname, restart)
{ }


void MultiBunchDump::fillHeader() {

    if (this->hasColumns()) {
        return;
    }

    columns_m.addColumn("t", "double", "ns", "Time");
    columns_m.addColumn("s", "double", "m", "Path length");
    columns_m.addColumn("azimuth", "double", "deg", "Azimuth in global coordinates");
    columns_m.addColumn("radius", "double", "mm", "Radius in global coordinates");
    columns_m.addColumn("numParticles", "long", "1", "Number of Macro Particles");
    columns_m.addColumn("energy", "double", "MeV", "Mean Bunch Energy");
    columns_m.addColumn("dE", "double", "MeV", "energy spread of the beam");

    columns_m.addColumn("rms_x", "double", "m", "RMS Beamsize in x");
    columns_m.addColumn("rms_y", "double", "m", "RMS Beamsize in y");
    columns_m.addColumn("rms_s", "double", "m", "RMS Beamsize in s");

    columns_m.addColumn("rms_px", "double", "1", "RMS Normalized Momenta in x");
    columns_m.addColumn("rms_py", "double", "1", "RMS Normalized Momenta in y");
    columns_m.addColumn("rms_ps", "double", "1", "RMS Normalized Momenta in s");

    columns_m.addColumn("emit_x", "double", "m", "Normalized Emittance x");
    columns_m.addColumn("emit_y", "double", "m", "Normalized Emittance y");
    columns_m.addColumn("emit_s", "double", "m", "Normalized Emittance s");

    columns_m.addColumn("mean_x", "double", "m", "Mean Beam Position in x");
    columns_m.addColumn("mean_y", "double", "m", "Mean Beam Position in y");
    columns_m.addColumn("mean_s", "double", "m", "Mean Beam Position in s");

    columns_m.addColumn("xpx", "double", "1", "RMS Correlation xpx");
    columns_m.addColumn("ypy", "double", "1", "RMS Correlation ypy");
    columns_m.addColumn("zpz", "double", "1", "RMS Correlation zpz");

    columns_m.addColumn("halo_x", "double", "1", "Halo in x");
    columns_m.addColumn("halo_y", "double", "1", "Halo in y");
    columns_m.addColumn("halo_z", "double", "1", "Halo in z");

    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;
    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Multi Bunch Statistics data '"
       << OpalData::getInstance()->getInputFn()
       << "' " << dateStr << " " << timeStr;

    this->addDescription(ss.str(), "multi bunch stat parameters");

    this->addDefaultParameters();

    this->addInfo("ascii", 1);
}


void MultiBunchDump::write(PartBunchBase<double, 3>* beam,
                           const beaminfo_t& binfo) {

    if ( Ippl::myNode() > 0)
        return;

    fillHeader();

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", binfo.time);
    columns_m.addColumnValue("s", binfo.pathlength);
    columns_m.addColumnValue("azimuth", binfo.azimuth);
    columns_m.addColumnValue("radius", binfo.radius);
    columns_m.addColumnValue("numParticles", binfo.nParticles);
    columns_m.addColumnValue("energy", binfo.ekin);
    columns_m.addColumnValue("dE", binfo.dEkin);

    columns_m.addColumnValue("rms_x", binfo.rrms[0]);
    columns_m.addColumnValue("rms_y", binfo.rrms[1]);
    columns_m.addColumnValue("rms_s", binfo.rrms[2]);

    columns_m.addColumnValue("rms_px", binfo.prms[0]);
    columns_m.addColumnValue("rms_py", binfo.prms[1]);
    columns_m.addColumnValue("rms_ps", binfo.prms[2]);

    columns_m.addColumnValue("emit_x", binfo.emit[0]);
    columns_m.addColumnValue("emit_y", binfo.emit[1]);
    columns_m.addColumnValue("emit_s", binfo.emit[2]);

    columns_m.addColumnValue("mean_x", binfo.mean[0]);
    columns_m.addColumnValue("mean_y", binfo.mean[1]);
    columns_m.addColumnValue("mean_s", binfo.mean[2]);

    columns_m.addColumnValue("xpx", binfo.correlation[0]);
    columns_m.addColumnValue("ypy", binfo.correlation[1]);
    columns_m.addColumnValue("zpz", binfo.correlation[2]);

    columns_m.addColumnValue("halo_x", binfo.halo[0]);
    columns_m.addColumnValue("halo_y", binfo.halo[1]);
    columns_m.addColumnValue("halo_z", binfo.halo[2]);

    this->writeRow();

    this->close();
}