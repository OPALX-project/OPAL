#include "MultiBunchDump.h"

#include <sstream>

#include "Utilities/Timer.h"
#include "Utilities/Util.h"

#include "Ippl.h"

extern Inform *gmsg;

MultiBunchDump::MultiBunchDump(const std::string& fname, bool restart,
                               const short& bunch)
    : SDDSWriter(fname, restart)
    , isFirst_m(true)
    , bunch_m(bunch)
{ }


void MultiBunchDump::fillHeader() {

    if ( !isFirst_m ) {
        return;
    }
    isFirst_m = false;
    
    columns_m.addColumn("t", "double", "ns", "Time");
    columns_m.addColumn("azimuth", "double", "deg", "Azimuth in global coordinates");
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
                           const double& azimuth) {
    // update binfo_m
    bool isOk = calcBunchBeamParameters(beam);
    binfo_m.azimuth = azimuth;

    if ( Ippl::myNode() > 0 || !isOk )
        return;

    fillHeader();

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", binfo_m.time);
    columns_m.addColumnValue("azimuth", binfo_m.azimuth);
    columns_m.addColumnValue("numParticles", binfo_m.nParticles);
    columns_m.addColumnValue("energy", binfo_m.ekin);
    columns_m.addColumnValue("dE", binfo_m.dEkin);

    columns_m.addColumnValue("rms_x", binfo_m.rrms[0]);
    columns_m.addColumnValue("rms_y", binfo_m.rrms[1]);
    columns_m.addColumnValue("rms_s", binfo_m.rrms[2]);

    columns_m.addColumnValue("rms_px", binfo_m.prms[0]);
    columns_m.addColumnValue("rms_py", binfo_m.prms[1]);
    columns_m.addColumnValue("rms_ps", binfo_m.prms[2]);

    columns_m.addColumnValue("emit_x", binfo_m.emit[0]);
    columns_m.addColumnValue("emit_y", binfo_m.emit[1]);
    columns_m.addColumnValue("emit_s", binfo_m.emit[2]);

    columns_m.addColumnValue("mean_x", binfo_m.mean[0]);
    columns_m.addColumnValue("mean_y", binfo_m.mean[1]);
    columns_m.addColumnValue("mean_s", binfo_m.mean[2]);

    columns_m.addColumnValue("halo_x", binfo_m.halo[0]);
    columns_m.addColumnValue("halo_y", binfo_m.halo[1]);
    columns_m.addColumnValue("halo_z", binfo_m.halo[2]);

    this->writeRow();

    this->close();
}


bool MultiBunchDump::calcBunchBeamParameters(PartBunchBase<double, 3>* beam) {
    if ( !OpalData::getInstance()->isInOPALCyclMode() ) {
        return false;
    }

    const unsigned long localNum = beam->getLocalNum();

    long int bunchTotalNum = 0;
    unsigned long bunchLocalNum = 0;

    /* container:
     *
     * ekin, <x>, <y>, <z>, <p_x>, <p_y>, <p_z>,
     * <x^2>, <y^2>, <z^2>, <p_x^2>, <p_y^2>, <p_z^2>,
     * <xp_x>, <y_py>, <zp_z>,
     * <x^3>, <y^3>, <z^3>, <x^4>, <y^4>, <z^4>
     */
    const unsigned int dim = PartBunchBase<double, 3>::Dimension;
    std::vector<double> local(7 * dim + 1);

    for(unsigned long k = 0; k < localNum; ++ k) {
        if ( beam->bunchNum[k] != bunch_m ) { //|| ID[k] == 0 ) {
            continue;
        }

        ++bunchTotalNum;
        ++bunchLocalNum;

        // ekin
        local[0] += std::sqrt(dot(beam->P[k], beam->P[k]) + 1.0);

        for (unsigned int i = 0; i < dim; ++i) {

            double r = beam->R[k](i);
            double p = beam->P[k](i);

            // <x>, <y>, <z>
            local[i + 1] += r;

            // <p_x>, <p_y, <p_z>
            local[i + dim + 1] += p;

            // <x^2>, <y^2>, <z^2>
            double r2 = r * r;
            local[i + 2 * dim + 1] += r2;

            // <p_x^2>, <p_y^2>, <p_z^2>
            local[i + 3 * dim + 1] += p * p;

            // <xp_x>, <y_py>, <zp_z>
            local[i + 4 * dim + 1] += r * p;

            // <x^3>, <y^3>, <z^3>
            local[i + 5 * dim + 1] += r2 * r;

            // <x^4>, <y^4>, <z^4>
            local[i + 6 * dim + 1] += r2 * r2;
        }
    }

    // inefficient
    allreduce(bunchTotalNum, 1, std::plus<long int>());

    if ( bunchTotalNum == 0 )
        return false;

    // ekin
    const double m0 = beam->getM() * 1.0e-6;
    local[0] -= bunchLocalNum;
    local[0] *= m0;

    allreduce(local.data(), local.size(), std::plus<double>());

    double invN = 1.0 / double(bunchTotalNum);
    binfo_m.ekin = local[0] * invN;

    binfo_m.time       = beam->getT() * 1e9;  // ns
    binfo_m.nParticles = bunchTotalNum;

    for (unsigned int i = 0; i < dim; ++i) {

        double w = local[i + 1] * invN;
        double pw = local[i + dim + 1] * invN;
        double w2 = local[i + 2 * dim + 1] * invN;
        double pw2 = local[i + 3 * dim + 1] * invN;
        double wpw = local[i + 4 * dim + 1] * invN;
        double w3 = local[i + 5 * dim + 1] * invN;
        double w4 = local[i + 6 * dim + 1] * invN;

        // <x>, <y>, <z>
        binfo_m.mean[i] = w;

        // central: <p_w^2> - <p_w>^2 (w = x, y, z)
        binfo_m.prms[i] = pw2 - pw * pw;
        if ( binfo_m.prms[i] < 0 ) {
            binfo_m.prms[i] = 0.0;
        }

        // central: <wp_w> - <w><p_w>
        wpw = wpw - w * pw;

        // central: <w^2> - <w>^2 (w = x, y, z)
        binfo_m.rrms[i] = w2 - w * w;

        // central: normalized emittance
        binfo_m.emit[i] = (binfo_m.rrms[i] * binfo_m.prms[i] - wpw * wpw);
        binfo_m.emit[i] =  std::sqrt(std::max(binfo_m.emit[i], 0.0));

        // central: <w^4> - 4 * <w> * <w^3> + 6 * <w>^2 * <w^2> - 3 * <w>^4
        double tmp = w4
                   - 4.0 * w * w3
                   + 6.0 * w * w * w2
                   - 3.0 * w * w * w * w;
        binfo_m.halo[i] = tmp / ( binfo_m.rrms[i] * binfo_m.rrms[i] );

        // central: sqrt(<w^2> - <w>^2) (w = x, y, z)
        binfo_m.rrms[i] = std::sqrt(binfo_m.rrms[i]);

        // central: sqrt(<p_w^2> - <p_w>^2)
        binfo_m.prms[i] = std::sqrt(binfo_m.prms[i]);
    }

    double tmp = 1.0 / std::pow(binfo_m.ekin / m0 + 1.0, 2.0);
    binfo_m.dEkin = binfo_m.prms[1] * m0 * std::sqrt(1.0 - tmp);

    return true;
}
