#include "MultiBunchHandler.h"

#include "Structure/H5PartWrapperForPC.h"

//FIXME Remove headers and dynamic_cast in
#include "Algorithms/PartBunch.h"
#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

extern Inform *gmsg;

MultiBunchHandler::MultiBunchHandler(PartBunchBase<double, 3> *beam,
                                     const int& numBunch,
                                     const double& eta,
                                     const double& para,
                                     const std::string& mode,
                                     const std::string& binning)
    : onebunch_m(OpalData::getInstance()->getInputBasename() + "-onebunch.h5")
    , numBunch_m(numBunch)
    , eta_m(eta)
    , coeffDBunches_m(para)
    , radiusLastTurn_m(0.0)
    , radiusThisTurn_m(0.0)
    , bunchCount_m(1)
{
    PAssert_GT(numBunch, 1);

    binfo_m.reserve(numBunch);
    for (int i = 0; i < beam->getNumBunch(); ++i) {
        binfo_m.push_back(beaminfo_t());
    }

    this->setBinning(binning);

    // mode of generating new bunches:
    // "FORCE" means generating one bunch after each revolution, until get "TURNS" bunches.
    // "AUTO" means only when the distance between two neighbor bunches is below the limitation,
    //        then starts to generate new bunches after each revolution,until get "TURNS" bunches;
    //        otherwise, run single bunch track

    *gmsg << "***---------------------------- MULTI-BUNCHES MULTI-ENERGY-BINS MODE "
            << "----------------------------*** " << endl;

    // only for regular  run of multi bunches, instantiate the  PartBins class
    // note that for restart run of multi bunches, PartBins class is instantiated in function
    // Distribution::doRestartOpalCycl()
    if (!OpalData::getInstance()->inRestartRun()) {

        // already exist bins number initially
        const int BinCount = 1;
        //allowed maximal bin
        const int MaxBinNum = 1000;

        // initialize particles number for each bin (both existed and not yet emmitted)
        size_t partInBin[MaxBinNum] = {0};
        partInBin[0] =  beam->getTotalNum();

        beam->setPBins(new PartBinsCyc(MaxBinNum, BinCount, partInBin));
        // the allowed maximal bin number is set to 100
        //beam->setPBins(new PartBins(100));

        this->setMode(mode);

    } else {
        if(beam->pbin_m->getLastemittedBin() < 2) {
            *gmsg << "In this restart job, the multi-bunches mode is forcely set to AUTO mode." << endl;
            mode_m = MB_MODE::AUTO;
        } else {
            *gmsg << "In this restart job, the multi-bunches mode is forcely set to FORCE mode." << endl
                    << "If the existing bunch number is less than the specified number of TURN, "
                    << "readin the phase space of STEP#0 from h5 file consecutively" << endl;
            mode_m = MB_MODE::FORCE;
        }
    }
}


void MultiBunchHandler::saveBunch(PartBunchBase<double, 3> *beam)
{
    static IpplTimings::TimerRef saveBunchTimer = IpplTimings::getTimer("Save Bunch H5");
    IpplTimings::startTimer(saveBunchTimer);
    *gmsg << endl;
    *gmsg << "* Store beam to H5 file for multibunch simulation ... ";

    Ppos_t coord, momentum;
    ParticleAttrib<double> mass, charge;
    ParticleAttrib<short> ptype;

    std::size_t localNum = beam->getLocalNum();

    coord.create(localNum);
    coord = beam->R;

    momentum.create(localNum);
    momentum = beam->P;

    mass.create(localNum);
    mass = beam->M;

    charge.create(localNum);
    charge = beam->Q;

    ptype.create(localNum);
    ptype = beam->PType;

    std::map<std::string, double> additionalAttributes = {
        std::make_pair("REFPR", 0.0),
        std::make_pair("REFPT", 0.0),
        std::make_pair("REFPZ", 0.0),
        std::make_pair("REFR", 0.0),
        std::make_pair("REFTHETA", 0.0),
        std::make_pair("REFZ", 0.0),
        std::make_pair("AZIMUTH", 0.0),
        std::make_pair("ELEVATION", 0.0),
        std::make_pair("B-ref_x",  0.0),
        std::make_pair("B-ref_z",  0.0),
        std::make_pair("B-ref_y",  0.0),
        std::make_pair("E-ref_x",  0.0),
        std::make_pair("E-ref_z",  0.0),
        std::make_pair("E-ref_y",  0.0),
        std::make_pair("B-head_x", 0.0),
        std::make_pair("B-head_z", 0.0),
        std::make_pair("B-head_y", 0.0),
        std::make_pair("E-head_x", 0.0),
        std::make_pair("E-head_z", 0.0),
        std::make_pair("E-head_y", 0.0),
        std::make_pair("B-tail_x", 0.0),
        std::make_pair("B-tail_z", 0.0),
        std::make_pair("B-tail_y", 0.0),
        std::make_pair("E-tail_x", 0.0),
        std::make_pair("E-tail_z", 0.0),
        std::make_pair("E-tail_y", 0.0)
    };

    H5PartWrapperForPC h5wrapper(onebunch_m, H5_O_WRONLY);
    h5wrapper.writeHeader();
    h5wrapper.writeStep(beam, additionalAttributes);
    h5wrapper.close();

    *gmsg << "Done." << endl;
    IpplTimings::stopTimer(saveBunchTimer);
}


bool MultiBunchHandler::readBunch(PartBunchBase<double, 3> *beam,
                                  const PartData& ref)
{
    static IpplTimings::TimerRef readBunchTimer = IpplTimings::getTimer("Read Bunch H5");
    IpplTimings::startTimer(readBunchTimer);
    *gmsg << endl;
    *gmsg << "* Read beam from H5 file for multibunch simulation ... ";

    std::size_t localNum = beam->getLocalNum();

    /*
     * 2nd argument: 0  --> step zero is fine since the file has only this step
     * 3rd argument: "" --> onebunch_m is used
     * 4th argument: H5_O_RDONLY does not work with this constructor
     */
    H5PartWrapperForPC h5wrapper(onebunch_m, 0, "", H5_O_WRONLY);

    size_t numParticles = h5wrapper.getNumParticles();

    const int bunchNum = bunchCount_m - 1;

    beam->setTotalNumPerBunch(numParticles, bunchNum);

    if ( numParticles == 0 ) {
        throw OpalException("MultiBunchHandler::readBunch()",
                            "No particles in file " + onebunch_m + ".");
    }

    size_t numParticlesPerNode = numParticles / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;

    PAssert_LT(firstParticle, lastParticle +1);

    numParticles = lastParticle - firstParticle + 1;

    beam->setLocalNumPerBunch(numParticles, bunchNum);

    //FIXME
    std::unique_ptr<PartBunchBase<double, 3> > tmpBunch = nullptr;
#ifdef ENABLE_AMR
    AmrPartBunch* amrbunch_p = dynamic_cast<AmrPartBunch*>(beam);
    if ( amrbunch_p != nullptr ) {
        tmpBunch.reset(new AmrPartBunch(&ref,
                                        amrbunch_p->getAmrParticleBase()));
    } else
#endif
        tmpBunch.reset(new PartBunch(&ref));

    tmpBunch->create(numParticles);

    h5wrapper.readStep(tmpBunch.get(), firstParticle, lastParticle);
    h5wrapper.close();

    beam->create(numParticles);

    for(size_t ii = 0; ii < numParticles; ++ ii, ++ localNum) {
        beam->R[localNum] = tmpBunch->R[ii];
        beam->P[localNum] = tmpBunch->P[ii];
        beam->M[localNum] = tmpBunch->M[ii];
        beam->Q[localNum] = tmpBunch->Q[ii];
        beam->PType[localNum] = ParticleType::REGULAR;
        beam->Bin[localNum] = bunchNum;
        beam->bunchNum[localNum] = bunchNum;
    }

    beam->boundp();

    binfo_m.push_back(beaminfo_t(injection_m));

    *gmsg << "Done." << endl;

    IpplTimings::stopTimer(readBunchTimer);
    return true;
}


short MultiBunchHandler::injectBunch(PartBunchBase<double, 3> *beam,
                                     const PartData& ref,
                                     bool& flagTransition)
{
    short result = 0;
    if ((bunchCount_m == 1) && (mode_m == MB_MODE::AUTO) && (!flagTransition)) {

        // we have still a single bunch
        beam->setTotalNumPerBunch(beam->getTotalNum(), 0);
        beam->setLocalNumPerBunch(beam->getLocalNum(), 0);

        // If all of the following conditions are met, this code will be executed
        // to check the distance between two neighboring bunches:
        // 1. Only one bunch exists (bunchCount_m == 1)
        // 2. We are in multi-bunch mode, AUTO sub-mode (mode_m == 2)
        // 3. It has been a full revolution since the last check (stepsNextCheck)

        *gmsg << "* MBM: Checking for automatically injecting new bunch ..." << endl;

        //beam->R *= Vector_t(0.001); // mm --> m
        beam->calcBeamParameters();
        //beam->R *= Vector_t(1000.0); // m --> mm

        Vector_t Rmean = beam->get_centroid(); // m

        radiusThisTurn_m = std::hypot(Rmean[0],Rmean[1]);

        Vector_t Rrms = beam->get_rrms(); // m

        double XYrms = std::hypot(Rrms[0], Rrms[1]);

        // If the distance between two neighboring bunches is less than 5 times of its 2D rms size
        // start multi-bunch simulation, fill current phase space to initialR and initialP arrays
        if ((radiusThisTurn_m - radiusLastTurn_m) < coeffDBunches_m * XYrms) {
            // since next turn, start multi-bunches
            saveBunch(beam);
            flagTransition = true;
        }

        *gmsg << "* MBM: RLastTurn = " << radiusLastTurn_m << " [m]" << endl;
        *gmsg << "* MBM: RThisTurn = " << radiusThisTurn_m << " [m]" << endl;
        *gmsg << "* MBM: XYrms = " << XYrms    << " [m]" << endl;

        radiusLastTurn_m = radiusThisTurn_m;
        result = 1;
    }

    else if (bunchCount_m < numBunch_m) {
        // Matthias: SteptoLastInj was used in MtsTracker, removed by DW in GenericTracker

        // If all of the following conditions are met, this code will be executed
        // to read new bunch from hdf5 format file:
        // 1. We are in multi-bunch mode (numBunch_m > 1)
        // 2. It has been a full revolution since the last check
        // 3. Number of existing bunches is less than the desired number of bunches
        // 4. FORCE mode, or AUTO mode with flagTransition = true
        // Note: restart from 1 < BunchCount < numBunch_m must be avoided.
        *gmsg << "* MBM: Injecting a new bunch ..." << endl;

        bunchCount_m++;

        beam->setNumBunch(bunchCount_m);

        // read initial distribution from h5 file
        switch ( mode_m ) {
            case MB_MODE::FORCE:
            case MB_MODE::AUTO:
                readBunch(beam, ref);
                updateParticleBins(beam);
                calcBunchBeamParameters(beam, bunchCount_m - 1);
                break;
            default:
                throw OpalException("MultiBunchHandler::injectBunch()",
                                    "We shouldn't be here in single bunch mode.");
        }

        Ippl::Comm->barrier();

        *gmsg << "* MBM: Bunch " << bunchCount_m
              << " injected, total particle number = "
              << beam->getTotalNum() << endl;
        result = 2;
    }
    return result;
}


void MultiBunchHandler::updateParticleBins(PartBunchBase<double, 3> *beam) {
    if (bunchCount_m < 2)
        return;

    static IpplTimings::TimerRef binningTimer = IpplTimings::getTimer("Particle Binning");
    IpplTimings::startTimer(binningTimer);
    switch ( binning_m ) {
        case MB_BINNING::GAMMA:
            beam->resetPartBinID2(eta_m);
            break;
        case MB_BINNING::BUNCH:
            beam->resetPartBinBunch();
            break;
        default:
            beam->resetPartBinID2(eta_m);
    }
    IpplTimings::stopTimer(binningTimer);
}


void MultiBunchHandler::setMode(const std::string& mbmode) {
    if ( mbmode.compare("FORCE") == 0 ) {
        *gmsg << "FORCE mode: The multi bunches will be injected consecutively" << endl
              << "            after each revolution, until get \"TURNS\" bunches." << endl;
        mode_m = MB_MODE::FORCE;
    } else if ( mbmode.compare("AUTO") == 0 ) {
        *gmsg << "AUTO mode: The multi bunches will be injected only when the" << endl
              << "           distance between two neighboring bunches is below" << endl
              << "           the limitation. The control parameter is set to "
              << coeffDBunches_m << endl;
        mode_m = MB_MODE::AUTO;
    } else
        throw OpalException("MultiBunchHandler::setMode()",
                            "MBMODE name \"" + mbmode + "\" unknown.");
}


void MultiBunchHandler::setBinning(std::string binning) {

    binning = Util::toUpper(binning);

    if ( binning.compare("BUNCH") == 0 ) {
        *gmsg << "Use 'BUNCH' injection for binnning." << endl;
        binning_m = MB_BINNING::BUNCH;
    } else if ( binning.compare("GAMMA") == 0 ) {
        *gmsg << "Use 'GAMMA' for binning." << endl;
        binning_m = MB_BINNING::GAMMA;
    } else {
        throw OpalException("MultiBunchHandler::setBinning()",
                            "MB_BINNING name \"" + binning + "\" unknown.");
    }
}


void MultiBunchHandler::setRadiusTurns(const double& radius) {
    if ( mode_m != MB_MODE::AUTO )
        return;

    radiusLastTurn_m = radius;
    radiusThisTurn_m = radiusLastTurn_m;

    if (OpalData::getInstance()->inRestartRun()) {
        *gmsg << "Radial position at restart position = ";
    } else {
        *gmsg << "Initial radial position = ";
    }
    // New OPAL 2.0: Init in m -DW
    *gmsg << radiusThisTurn_m << " m" << endl;
}


bool MultiBunchHandler::calcBunchBeamParameters(PartBunchBase<double, 3>* beam,
                                                short bunchNr)
{
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

    beaminfo_t& binfo = getBunchInfo(bunchNr);

    for(unsigned long k = 0; k < localNum; ++k) {
        if ( beam->bunchNum[k] != bunchNr ) { //|| ID[k] == 0 ) {
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

    // here we also update the number of particles of *this* bunch
    if (bunchNr >= (short)beam->getNumBunch())
        throw OpalException("MultiBunchHandler::calcBunchBeamParameters()",
                            "Bunch number " + std::to_string(bunchNr) +
                            " exceeds bunch index " + std::to_string(beam->getNumBunch() - 1));

    beam->setTotalNumPerBunch(bunchTotalNum, bunchNr);
    beam->setLocalNumPerBunch(bunchLocalNum, bunchNr);

    if ( bunchTotalNum == 0 )
        return false;

    // ekin
    const double m0 = beam->getM() * 1.0e-6;
    local[0] -= bunchLocalNum;
    local[0] *= m0;

    allreduce(local.data(), local.size(), std::plus<double>());

    double invN = 1.0 / double(bunchTotalNum);
    binfo.ekin = local[0] * invN;

    binfo.time       = beam->getT() * 1e9;  // ns
    binfo.nParticles = bunchTotalNum;

    for (unsigned int i = 0; i < dim; ++i) {

        double w   = local[i + 1] * invN;
        double pw  = local[i + dim + 1] * invN;
        double w2  = local[i + 2 * dim + 1] * invN;
        double pw2 = local[i + 3 * dim + 1] * invN;
        double wpw = local[i + 4 * dim + 1] * invN;
        double w3  = local[i + 5 * dim + 1] * invN;
        double w4  = local[i + 6 * dim + 1] * invN;

        // <x>, <y>, <z>
        binfo.mean[i] = w;

        // central: <p_w^2> - <p_w>^2 (w = x, y, z)
        binfo.prms[i] = pw2 - pw * pw;
        if ( binfo.prms[i] < 0 ) {
            binfo.prms[i] = 0.0;
        }

        // central: <wp_w> - <w><p_w>
        wpw = wpw - w * pw;

        // central: <w^2> - <w>^2 (w = x, y, z)
        binfo.rrms[i] = w2 - w * w;

        // central: normalized emittance
        binfo.emit[i] = (binfo.rrms[i] * binfo.prms[i] - wpw * wpw);
        binfo.emit[i] =  std::sqrt(std::max(binfo.emit[i], 0.0));

        // central: <w^4> - 4 * <w> * <w^3> + 6 * <w>^2 * <w^2> - 3 * <w>^4
        double tmp = w4
                   - 4.0 * w * w3
                   + 6.0 * w * w * w2
                   - 3.0 * w * w * w * w;
        binfo.halo[i] = tmp / ( binfo.rrms[i] * binfo.rrms[i] );

        // central: sqrt(<w^2> - <w>^2) (w = x, y, z)
        binfo.rrms[i] = std::sqrt(binfo.rrms[i]);

        // central: sqrt(<p_w^2> - <p_w>^2)
        binfo.prms[i] = std::sqrt(binfo.prms[i]);

        // central: rms correlation --> (<wp_w> - <w><p_w>) / sqrt(<w^2> * <p_w^2>)
        double denom = 1.0 / (binfo.rrms[i] * binfo.prms[i]);
        binfo.correlation[i] = (std::isfinite(denom)) ? wpw * denom : 0.0;
    }

    double tmp = 1.0 / std::pow(binfo.ekin / m0 + 1.0, 2.0);
    binfo.dEkin = binfo.prms[1] * m0 * std::sqrt(1.0 - tmp);

    return true;
}


void MultiBunchHandler::updateTime(const double& dt) {
    for (short b = 0; b < bunchCount_m; ++b) {
        binfo_m[b].time += dt;
    }
}


void MultiBunchHandler::updatePathLength(const std::vector<double>& lpaths) {
    PAssert_EQ(bunchCount_m, (short)lpaths.size() - 1);
    for (short b = 0; b < bunchCount_m; ++b) {
        binfo_m[b].pathlength += lpaths[b];
    }
}
