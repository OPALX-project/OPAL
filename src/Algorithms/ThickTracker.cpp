// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickTracker
//   The visitor class for building a map of given order for a beamline
//   using a finite-length lenses for all elements.
//   Multipole-like elements are done by expanding the Lie series.
//
// ------------------------------------------------------------------------
//
// $Author: ganz_p $
//
// ------------------------------------------------------------------------


#include <cfloat>
#include <fstream>
#include <typeinfo>

#include "Algorithms/ThickTracker.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"

#include "Classic/Algorithms/PartData.h"

#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"

#include "Physics/Physics.h"

//
// Class ThickTracker
// ------------------------------------------------------------------------
//
ThickTracker::ThickTracker(const Beamline &beamline,
                           const PartData &reference,
                           bool revBeam, bool revTrack)
    : Tracker(beamline, reference, revBeam, revTrack)
    , hamiltonian_m(1)
    , RefPartR_m(0.0)
    , RefPartP_m(0.0)
    , itsDataSink_m(nullptr)
    , itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection())
    , zstart_m(0.0)
    , zstop_m(0.0)
    , threshold_m(1.0e-6)
    , truncOrder_m(1) // linear
    , mapCreation_m(   IpplTimings::getTimer("mapCreation"))
    , mapCombination_m(IpplTimings::getTimer("mapCombination"))
    , mapTracking_m(   IpplTimings::getTimer("mapTracking"))
{
    CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
                                   beamline.getInitialDirection());
    referenceToLabCSTrafo_m = labToRef.inverted();
}


ThickTracker::ThickTracker(const Beamline &beamline,
                           PartBunchBase<double, 3> *bunch,
                           Beam &beam,
                           DataSink &ds,
                           const PartData &reference,
                           bool revBeam, bool revTrack,
                           const std::vector<unsigned long long> &maxSteps,
                           double zstart,
                           const std::vector<double> &zstop,
                           const std::vector<double> &dt,
                           const int& truncOrder)
    : Tracker(beamline, bunch, reference, revBeam, revTrack)
    , hamiltonian_m(truncOrder)
    , RefPartR_m(0.0)
    , RefPartP_m(0.0)
    , itsDataSink_m(&ds)
    , itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection())
    , zstart_m(zstart)
    , zstop_m(zstop[0])
    , threshold_m(1.0e-6)
    , truncOrder_m(truncOrder)
    , mapCreation_m(   IpplTimings::getTimer("mapCreation"))
    , mapCombination_m(IpplTimings::getTimer("mapCombination"))
    , mapTracking_m(   IpplTimings::getTimer("mapTracking"))
{
    if ( zstop.size() > 1 )
        throw OpalException("ThickTracker::ThickTracker()",
                            "Multiple tracks not yet supported.");


    CoordinateSystemTrafo labToRef(beamline.getOrigin3D(),
                                   beamline.getInitialDirection());
    referenceToLabCSTrafo_m = labToRef.inverted();
}


ThickTracker::~ThickTracker()
{}



void ThickTracker::visitBeamline(const Beamline &bl) {

    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    if (fbl->getRelativeFlag()) {
        *gmsg << " do stuff" << endl;
        OpalBeamline stash(fbl->getOrigin3D(), fbl->getInitialDirection());
        stash.swap(itsOpalBeamline_m);
        fbl->iterate(*this, false);
        itsOpalBeamline_m.prepareSections();
        itsOpalBeamline_m.compute3DLattice();
        stash.merge(itsOpalBeamline_m);
        stash.swap(itsOpalBeamline_m);
    } else {
        fbl->iterate(*this, false);
    }
}


void ThickTracker::prepareSections() {
    itsBeamline_m.accept(*this);
    itsOpalBeamline_m.prepareSections();
}




/*
//TODO complete and test fringe fields
void ThickTracker::insertFringeField(SBend* pSBend, lstruct_t& mBL,
        double& beta0, double& gamma0, double& P0, double& q, std::array<double,2>& paramFringe, std::string e){


    series_t H;
    map_t tempMap, fieldMap;
    lstruct_t::iterator mBLit;
    double lenFringe = paramFringe[0] + paramFringe[1];

    //double intHeight = 0.1; //:TODO: change or argument for that value
    double stepSize = 0.01; //:TODO: change or argument for that value maybe minStepIfSC?


    //entrFringe = std::fabs(intHeight * std::tan(pSBend->getEntranceAngle()));

    //get proper slices for the fringe field
    double nSlices = std::ceil(2 * lenFringe / stepSize);
    stepSize = 2 * lenFringe / nSlices;

    //Create structMapTracking and push in BeamLine
    double K0= pSBend ->getB()*(Physics::c/itsReference.getP());
    K0= std::round(K0*1e6)/1e6 *q*(Physics::c/P0);

    double h = 1./ pSBend ->getBendRadius();               //inverse bending radius [1/m]
    K0<0 ? h = -1*h : h = 1*h;

    series_t az = - K0 / (2 * lenFringe) * (hamiltonian_m.y*hamiltonian_m.y - hamiltonian_m.x*hamiltonian_m.x) * std::tan(pSBend->getEntranceAngle());
    if (e == "out") {
            az = -az;
        }
    std::cout << az.getMaxOrder() << std::endl;
    for ( int i = 0 ; i < nSlices; i++ ){
        double l = stepSize * i;
        series_t ax = 0.5 * K0 / (2 * lenFringe) * (l*l - hamiltonian_m.y*hamiltonian_m.y);
        H = hamiltonian_m.bendFringe(beta0, gamma0, h, K0, ax, az);
        tempMap = ExpMap(-H * stepSize ,truncOrder_m);
        fieldMap = (tempMap * fieldMap).truncate(truncOrder_m);
    }

    structMapTracking fringeField;
    fringeField.elementName=pSBend->getName() + "accumulated Fringe Field Map";
    if (e == "in"){
        fringeField.elementPos=(pSBend->getElementPosition() - paramFringe[0]);
    } else if (e == "out") {
        fringeField.elementPos=(pSBend->getElementPosition() + pSBend->getArcLength() - paramFringe[1]);
    }
    fringeField.stepSize= lenFringe;
    fringeField.nSlices=1;
    fringeField.elementMap=fieldMap;
    mBL.push_back(fringeField);
}
*/


/**
 * @brief Algorithm for Thick Map-Tracking
 *
 */
void ThickTracker::execute() {

    /*
     * global settings
     */
    Inform msg("ThickTracker", *gmsg);

    OpalData::getInstance()->setInPrepState(true);

    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    if ( OpalData::getInstance()->hasPriorTrack() ||
         OpalData::getInstance()->inRestartRun() )
    {
        OpalData::getInstance()->setOpenMode(OpalData::OPENMODE::APPEND);
    }

    prepareSections();

    msg << "Truncation order: " << this->truncOrder_m << endl;


    this->checkElementOrder_m();

    this->fillGaps_m();

    //FIXME in local system only
    RefPartR_m = Vector_t(0.0);
    RefPartP_m = euclidean_norm(itsBunch_m->get_pmean_Distribution()) * Vector_t(0, 0, 1);

    itsBunch_m->set_sPos(zstart_m);

//     IpplTimings::startTimer(mapCreation_m);

    // TODO need a way to implement Initial dispersion Values
    //dispInitialVal[0][0]=  0.5069938765;
    //dispInitialVal[0][1]= -0.1681363086;
    OpalData::getInstance()->setInPrepState(false);

    track_m();

    // fMatrix_t sigMatrix = itsBunch_m->getSigmaMatrix();
    // linSigAnalyze(sigMatrix);
}


void ThickTracker::checkElementOrder_m() {

    // check order of beam line
    FieldList elements = itsOpalBeamline_m.getElementByType(ElementBase::ANY);
    beamline_t::const_iterator el = elements_m.cbegin();

    double currentEnd = zstart_m;
    for (FieldList::iterator it = elements.begin(); it != elements.end(); ++it) {
        double pos = it->getElement()->getElementPosition();
        // we have to take this length due to dipole --> arclength
        if ( currentEnd - pos > threshold_m ) {

            throw OpalException("ThickTracker::checkOrder_m()",
                                std::string("Elements are not in ascending order or overlap.") +
                                " Element Name: " + it->getElement()->getName() +
                                " ELEMEDGE position: " + std::to_string(pos) +
                                " Beamline end " + std::to_string(currentEnd) +
                                " element length " + std::to_string(std::get<2>(*el))) ;
        }
        currentEnd = pos + std::get<2>(*el);
        ++el;
        while (std::get<2>(*el) == 0) ++el; // skip zero-length elements (aka fringe fields)
    }
}


///Fills undefined beam path with a Drift Space
/** \f[H_{Drift}= \frac{\delta}{\beta_0} -
* \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } \f]
*/
void ThickTracker::fillGaps_m() {

    beamline_t tmp;

    FieldList elements = itsOpalBeamline_m.getElementByType(ElementBase::ANY);
    beamline_t::const_iterator el = elements_m.cbegin();

    double currentEnd = zstart_m;
    for (FieldList::iterator it = elements.begin(); it != elements.end(); ++it) {

        tmp.push_back( *el );

        double pos = it->getElement()->getElementPosition();

        double length = std::abs(pos - currentEnd);
        if ( length  > threshold_m ) {
            //FIXME if gamma changes this is an error!!!
            double gamma = itsReference.getGamma();
            tmp.push_back(std::make_tuple(hamiltonian_m.drift(gamma), 1, length));
        }

        currentEnd = pos + std::get<2>(*el);
        ++el;
    }

    double length = std::abs(zstop_m - currentEnd);
    if ( length > threshold_m ) {
        //FIXME if gamma changes this is an error!!!
        double gamma = itsReference.getGamma();
        tmp.push_back(std::make_tuple(hamiltonian_m.drift(gamma), 1, length));
    }

    elements_m.swap(tmp);
}


void ThickTracker::track_m()
{
    IpplTimings::startTimer(mapTracking_m);

    fMatrix_t tFMatrix;
    for (int i=0; i<6; i++){
        tFMatrix[i][i]=1.;
    }

    FMatrix<double, 1, 4> dispInitialVal;
    for (int i=2; i<4; i++){
        dispInitialVal[0][i]=0;
    }

    map_t transferMap;

    fMatrix_t refSigma;
    refSigma = (itsBunch_m->getSigmaMatrix());

    double spos = zstart_m;


    this->advanceDispersion_m(tFMatrix, dispInitialVal, spos);


    std::size_t step = 0;

    this->dump_m();

    //(1) Loop Beamline
    for(beamline_t::const_iterator it = elements_m.cbegin(); it != elements_m.end(); ++it) {
        //(2) Loop Slices

        const series_t& H          = std::get<0>(*it);
        const std::size_t& nSlices = std::get<1>(*it);
        const double& length       = std::get<2>(*it);

        const double ds = length / double(nSlices);

        map_t map = ExpMap(- H * ds, truncOrder_m);

        // global truncation order is truncOrder_m + 1
        map = map.truncate(truncOrder_m);

        for (std::size_t slice = 0; slice < nSlices; ++slice) {

            this->advanceParticles_m(map);

            this->concatenateMaps_m(map, transferMap);

            tFMatrix= transferMap.linearTerms();
            this->advanceDispersion_m(tFMatrix, dispInitialVal, spos);

            spos += ds;
            ++step;

            this->update_m(spos, step);

            this->dump_m();

            refSigma = itsBunch_m->getSigmaMatrix();

            //printPhaseShift(refSigma ,mapBeamLineit->elementMap.linearTerms(), N);
        }
    }


    this->write_m(transferMap);

    IpplTimings::stopTimer(mapTracking_m);
}


ThickTracker::particle_t
ThickTracker::particleToVector_m(const Vector_t& R,
                                 const Vector_t& P) const
{
    particle_t particle;
    for (int d = 0; d < 3; ++d) {
        particle[2 * d] = R(d);
        particle[2 *d + 1] = P(d);
    }
    return particle;
}


void ThickTracker::vectorToParticle_m(const particle_t& particle,
                                      Vector_t& R,
                                      Vector_t& P) const
{
    for (int d = 0; d < 3; ++d) {
        R(d) = particle[2 * d];
        P(d) = particle[2 *d + 1];
    }
}


void ThickTracker::write_m(const map_t& map) {

    if ( Ippl::myNode() == 0 ) {

        static bool first = true;

        std::string fn = OpalData::getInstance()->getInputBasename() + ".map";

        std::ofstream out;

        if ( first ) {
            first = false;
            out.open(fn, std::ios::out);
        } else {
            out.open(fn, std::ios::app);
        }

        out << std::setprecision(16)
            << map
            << std::endl;

        out.close();
    }
}


void ThickTracker::advanceParticles_m(const map_t& map) {
    for (std::size_t ip = 0; ip < itsBunch_m->getLocalNum(); ++ip) {

        particle_t particle = this->particleToVector_m(itsBunch_m->R[ip],
                                                       itsBunch_m->P[ip]);

        this->updateParticle_m(particle, map);

        this->vectorToParticle_m(particle,
                                 itsBunch_m->R[ip],
                                 itsBunch_m->P[ip]);
    }

    // update reference particle
    particle_t particle = this->particleToVector_m(RefPartR_m,
                                                   RefPartP_m);

    this->updateParticle_m(particle, map);

    this->vectorToParticle_m(particle,
                             RefPartR_m,
                             RefPartP_m);
}


void ThickTracker::updateParticle_m(particle_t& particle,
                                    const map_t& map)
{
    double betagamma = itsBunch_m->getInitialBeta() * itsBunch_m->getInitialGamma();

    //Units
    double pParticle= std::sqrt( particle[1] * particle[1] +
                                 particle[3] * particle[3] +
                                 particle[5] * particle[5]); // [beta gamma]

    particle[1] /= betagamma;
    particle[3] /= betagamma;

    particle[5] = std::sqrt( 1.0 + pParticle * pParticle ) / betagamma
                    - 1.0 / itsBunch_m->getInitialBeta();

    //Apply element map
    particle = map * particle;

    //Units back
    particle[1] *= betagamma;
    particle[3] *= betagamma;

    double tempGamma = (particle[5] + 1.0 / itsBunch_m->getInitialBeta())
            * betagamma ;

    pParticle =  std::sqrt( tempGamma * tempGamma -1.0);

    particle[5] = std::sqrt(pParticle * pParticle -
                            particle[1] * particle[1] -
                            particle[3] * particle[3] );
}


void ThickTracker::concatenateMaps_m(const map_t& x, map_t& y) {
    IpplTimings::startTimer(mapCombination_m);
    y = x * y;
    y = y.truncate(truncOrder_m);
    IpplTimings::stopTimer(mapCombination_m);
}


///Writes Dispersion in X and Y plane
/** used formula:
 * \f[\begin{pmatrix}\eta_{x} \\ \eta_{p_x}\end{pmatrix}_{s_1} =
 *   \begin{pmatrix} R_{11} & R_{12} \\ R_{21} & R_{22} \end{pmatrix}
 *  \cdot
 *  \begin{pmatrix} \eta_{x} \\ \eta_{p_x} \end{pmatrix}_{s_0}
 *  +
 *  \begin{pmatrix} R_{16} \\ R_{26} \end{pmatrix}\f]
 */
void ThickTracker::advanceDispersion_m(fMatrix_t tempMatrix,
                                       FMatrix<double, 1, 4> initialVal,
                                       double pos)
{
    if ( Ippl::myNode() == 0 ) {
        static bool first = true;

        std::ofstream out;
        std::string fn = OpalData::getInstance()->getInputBasename() + ".dispersion";

        if ( first ) {
            first = false;
            out.open(fn, std::ios::out);
        } else {
            out.open(fn, std::ios::app);
        }

        out << std::setprecision(16);


        FMatrix<double, 2, 2> subx, suby;
        FMatrix<double, 2, 1> subdx, subdy;
        FMatrix<double, 2, 1> dxi, dyi, dx, dy;

        for (int i=0; i<2; i++){
            dxi[i][0]= initialVal[0][i]; //
            dyi[i][0]= initialVal[0][i+2];

            subdx[i][0]=tempMatrix[i][5]*itsBunch_m->getInitialBeta();
            subdy[i][0]=tempMatrix[i+2][5]*itsBunch_m->getInitialBeta();

            for (int j=0; j<2; j++){

                subx[i][j]=tempMatrix[i][j];
                suby[i][j]=tempMatrix[i+2][j+2];

            }
        }

        dx= subx*dxi + subdx;
        dy= suby*dyi + subdy;



        out <<pos<< "\t" << dx[0][0] << "\t" << dx[0][1] << "\t" << dy[0][0] << "\t" << dy[0][1] << std::endl;
    }
}


void ThickTracker::dump_m() {

    if ( itsBunch_m->getTotalNum() == 0 )
        return;

    Inform msg("ThickTracker", *gmsg);

    msg << *itsBunch_m << endl;

    const std::size_t step = itsBunch_m->getGlobalTrackStep();

    bool psDump   = ((step + 1) % Options::psDumpFreq == 0);
    bool statDump = ((step + 1) % Options::statDumpFreq == 0);

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.

    Vector_t FDext[2];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    if (psDump || statDump) {

        Vector_t externalE, externalB;

        Vector_t rmin, rmax;
        itsBunch_m->get_bounds(rmin, rmax);

        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(referenceToLabCSTrafo_m.transformTo(RefPartR_m),
                                     referenceToLabCSTrafo_m.rotateTo(RefPartP_m),
                                     itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(),
                                     externalE,
                                     externalB);
        FDext[0] = referenceToLabCSTrafo_m.rotateFrom(externalB);
        FDext[1] = referenceToLabCSTrafo_m.rotateFrom(externalE * 1e-6);
    }


    if ( psDump ) {
        itsDataSink_m->dumpH5(itsBunch_m, FDext);
    }

    if ( statDump ) {
        std::vector<std::pair<std::string, unsigned int> > collimatorLosses;
        itsDataSink_m->dumpSDDS(itsBunch_m, FDext, collimatorLosses);
    }
}


void ThickTracker::update_m(const double& spos,
                            const std::size_t& step)
{
    // update dt and t
    double ds = spos - itsBunch_m->get_sPos();
    double gamma = itsBunch_m->get_gamma();
    double beta  = std::sqrt(1.0 - 1.0 / (gamma * gamma) );
    double dt = ds / Physics::c / beta;
    itsBunch_m->setdT(dt);
    itsBunch_m->setT(itsBunch_m->getT() + dt);

    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch_m->dt[i] = dt;
    }


    itsBunch_m->set_sPos(spos);
    itsBunch_m->setGlobalTrackStep(step);
    itsBunch_m->calcBeamParameters();
    itsBunch_m->calcEMean();
}