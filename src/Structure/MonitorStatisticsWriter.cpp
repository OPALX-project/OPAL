#include "MonitorStatisticsWriter.h"

#include "Structure/LossDataSink.h"
#include "Ippl.h"

MonitorStatisticsWriter::MonitorStatisticsWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }


void MonitorStatisticsWriter::fillHeader() {

    if (this->hasColumns()) {
        return;
    }

    this->addDescription("Statistics data of monitors",
                          "stat parameters");
    this->addDefaultParameters();

    columns_m.addColumn("name", "string", "", "Monitor name");
    columns_m.addColumn("s", "double", "m", "Longitudinal Position");
    columns_m.addColumn("t", "double", "ns", "Passage Time Reference Particle");
    columns_m.addColumn("numParticles", "long", "1", "Number of Macro Particles");
    columns_m.addColumn("rms_x", "double", "m", "RMS Beamsize in x");
    columns_m.addColumn("rms_y", "double", "m", "RMS Beamsize in y");
    columns_m.addColumn("rms_s", "double", "m", "RMS Beamsize in s");
    columns_m.addColumn("rms_t", "double", "ns", "RMS Passage Time");
    columns_m.addColumn("rms_px", "double", "1", "RMS Momenta in x");
    columns_m.addColumn("rms_py", "double", "1", "RMS Momenta in y");
    columns_m.addColumn("rms_ps", "double", "1", "RMS Momenta in s");
    columns_m.addColumn("emit_x", "double", "m", "Normalized Emittance x");
    columns_m.addColumn("emit_y", "double", "m", "Normalized Emittance y");
    columns_m.addColumn("emit_s", "double", "m", "Normalized Emittance s");
    columns_m.addColumn("mean_x", "double", "m", "Mean Beam Position in x");
    columns_m.addColumn("mean_y", "double", "m", "Mean Beam Position in y");
    columns_m.addColumn("mean_s", "double", "m", "Mean Beam Position in s");
    columns_m.addColumn("mean_t", "double", "ns", "Mean Passage Time");
    columns_m.addColumn("ref_x", "double", "m", "x coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_y", "double", "m", "y coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_z", "double", "m", "z coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_px", "double", "1", "x momentum of reference particle in lab cs");
    columns_m.addColumn("ref_py", "double", "1", "y momentum of reference particle in lab cs");
    columns_m.addColumn("ref_pz", "double", "1", "z momentum of reference particle in lab cs");
    columns_m.addColumn("min_x", "double", "m", "Max Beamsize in x");
    columns_m.addColumn("min_y", "double", "m", "Max Beamsize in y");
    columns_m.addColumn("min_s", "double", "m", "Max Beamsize in s");
    columns_m.addColumn("max_x", "double", "m", "Max Beamsize in x");
    columns_m.addColumn("max_y", "double", "m", "Max Beamsize in y");
    columns_m.addColumn("max_s", "double", "m", "Max Beamsize in s");
    columns_m.addColumn("xpx", "double", "1", "Correlation xpx");
    columns_m.addColumn("ypy", "double", "1", "Correlation ypy");
    columns_m.addColumn("zpz", "double", "1", "Correlation zpz");

    this->addInfo("ascii", 1);
}


void MonitorStatisticsWriter::addRow(const SetStatistics &set) {

    if ( Ippl::myNode() != 0 )
        return;


    this->fillHeader();

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("name", set.element_m);
    columns_m.addColumnValue("s", set.spos_m);
    columns_m.addColumnValue("t", set.refTime_m);
    columns_m.addColumnValue("numParticles", set.nTotal_m);
    columns_m.addColumnValue("rms_x", set.rrms_m(0));
    columns_m.addColumnValue("rms_y", set.rrms_m(1));
    columns_m.addColumnValue("rms_s", set.rrms_m(2));
    columns_m.addColumnValue("rms_t", set.trms_m * 1e9);
    columns_m.addColumnValue("rms_px", set.prms_m(0));
    columns_m.addColumnValue("rms_py", set.prms_m(1));
    columns_m.addColumnValue("rms_ps", set.prms_m(2));
    columns_m.addColumnValue("emit_x", set.eps_norm_m(0));
    columns_m.addColumnValue("emit_y", set.eps_norm_m(1));
    columns_m.addColumnValue("emit_s", set.eps_norm_m(2));
    columns_m.addColumnValue("mean_x", set.rmean_m(0));
    columns_m.addColumnValue("mean_y", set.rmean_m(1));
    columns_m.addColumnValue("mean_s", set.rmean_m(2));
    columns_m.addColumnValue("mean_t", set.tmean_m * 1e9);
    columns_m.addColumnValue("ref_x", set.RefPartR_m(0));
    columns_m.addColumnValue("ref_y", set.RefPartR_m(1));
    columns_m.addColumnValue("ref_z", set.RefPartR_m(2));
    columns_m.addColumnValue("ref_px", set.RefPartP_m(0));
    columns_m.addColumnValue("ref_py", set.RefPartP_m(1));
    columns_m.addColumnValue("ref_pz", set.RefPartP_m(2));
    columns_m.addColumnValue("min_x", set.rmin_m[0]);
    columns_m.addColumnValue("min_y", set.rmin_m[1]);
    columns_m.addColumnValue("min_s", set.rmin_m[2]);
    columns_m.addColumnValue("max_x", set.rmax_m[0]);
    columns_m.addColumnValue("max_y", set.rmax_m[1]);
    columns_m.addColumnValue("max_s", set.rmax_m[2]);
    columns_m.addColumnValue("xpx", set.rprms_m(0));
    columns_m.addColumnValue("ypy", set.rprms_m(1));
    columns_m.addColumnValue("zpz", set.rprms_m(2));

    this->writeRow();

    this->close();
}