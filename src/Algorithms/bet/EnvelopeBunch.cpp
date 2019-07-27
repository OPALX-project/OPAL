#include <iostream>
#include <cfloat>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>
#include <algorithm>
#include <typeinfo>

//FIXME: replace with IPPL Vector_t
#include <vector>
//FIXME: remove
#include <mpi.h>

#include "Algorithms/bet/math/root.h"     // root finding routines
#include "Algorithms/bet/math/linfit.h"   // linear fitting routines
#include "Algorithms/bet/math/savgol.h"   // savgol smoothing routine
#include "Algorithms/bet/math/rk.h"       // Runge-Kutta Integration

#include "Algorithms/bet/EnvelopeBunch.h"
#include "Utilities/Util.h"
#include "Utility/IpplTimings.h"

#define USE_HOMDYN_SC_MODEL

extern Inform *gmsg;

/** for space charge
    ignore for too low energies
    beta = 0.05 -> 640 eV
    beta = 0.10 ->   3 keV
    beta = 0.20 ->  11 keV
    beta = 0.30 ->  25 keV
    beta = 0.40 ->  47 keV
    beta = 0.50 ->  90 keV
    beta = 0.60 -> 128 keV
    beta = 0.70 -> 205 keV
    beta = 0.75 -> 261 keV
    beta = 0.80 -> 341 keV
    beta = 0.85 -> 460 keV
*/
#define BETA_MIN1 0.30     // minimum beta-value for space-charge calculations: start

#ifndef USE_HOMDYN_SC_MODEL
#define BETA_MIN2 0.45     // minimum beta-value for space-charge calculations: full impact
#endif

// Hack allows odeint in rk.C to be called with a class member function
static EnvelopeBunch *thisBunch = NULL;  // pointer to access calling bunch
static void Gderivs(double t, double Y[], double dYdt[]) { thisBunch->derivs(t, Y, dYdt); }
//static double Gcur(double z) { return thisBunch->currentProfile_m->get(z, itype_lin); }
static double rootValue = 0.0;

// used in setLShape for Gaussian
static void erfRoot(double x, double *fn, double *df) {
    double v = erfc(fabs(x));
    double eps = 1.0e-05;

    *fn = v - rootValue;
    *df = (erfc(fabs(x) + eps) - v) / eps;
}


EnvelopeBunch::EnvelopeBunch(const PartData *ref):
    PartBunch(ref),
    numSlices_m(0),
    numMySlices_m(0) {

    calcITimer_m = IpplTimings::getTimer("calcI");
    spaceChargeTimer_m = IpplTimings::getTimer("spaceCharge");
    isValid_m = true;

}


EnvelopeBunch::EnvelopeBunch(const std::vector<OpalParticle> &rhs,
                             const PartData *ref):
    PartBunch(ref),
    numSlices_m(0),
    numMySlices_m(0)
{}


EnvelopeBunch::~EnvelopeBunch() {
}

void EnvelopeBunch::calcBeamParameters() {
    Inform msg("calcBeamParameters");
    IpplTimings::startTimer(statParamTimer_m);
    double ex, ey, nex, ney, b0, bRms, bMax, bMin, g0, dgdt, gInc;
    double RxRms = 0.0, RyRms = 0.0, Px = 0.0, PxRms = 0.0, PxMax = 0.0, PxMin = 0.0, Py = 0.0, PyRms = 0.0, PyMax = 0.0, PyMin = 0.0;
    double Pz, PzMax, PzMin, PzRms;
    double x0Rms, y0Rms, zRms, zMax, zMin, I0, IRms, IMax, IMin;

    runStats(sp_beta, &b0, &bMax, &bMin, &bRms, &nValid_m);
    runStats(sp_I, &I0, &IMax, &IMin, &IRms, &nValid_m);
    runStats(sp_z, &z0_m, &zMax, &zMin, &zRms, &nValid_m);
    runStats(sp_Pz, &Pz, &PzMax, &PzMin, &PzRms, &nValid_m);

    if(solver_m & sv_radial) {
        runStats(sp_Rx, &Rx_m, &RxMax_m, &RxMin_m, &RxRms, &nValid_m);
        runStats(sp_Ry, &Ry_m, &RyMax_m, &RyMin_m, &RyRms, &nValid_m);
        runStats(sp_Px, &Px, &PxMax, &PxMin, &PxRms, &nValid_m);
        runStats(sp_Py, &Py, &PyMax, &PyMin, &PyRms, &nValid_m);
        calcEmittance(&nex, &ney, &ex, &ey, &nValid_m);
    }

    if(solver_m & sv_offaxis) {
        runStats(sp_x0, &x0_m, &x0Max_m, &x0Min_m, &x0Rms, &nValid_m);
        runStats(sp_y0, &y0_m, &y0Max_m, &y0Min_m, &y0Rms, &nValid_m);
    }

    calcEnergyChirp(&g0, &dgdt, &gInc, &nValid_m);
    double Bfz = AvBField();
    double Efz = AvEField();
    double mc2e = 1.0e-6 * Physics::EMASS * Physics::c * Physics::c / Physics::q_e;

    E_m = mc2e * (g0 - 1.0);
    dEdt_m = 1.0e-12 * mc2e * dgdt;
    Einc_m = mc2e * gInc;
    tau_m =  zRms / Physics::c;
    I_m = IMax;
    Irms_m = Q_m * nValid_m * Physics::c / (zRms * sqrt(Physics::two_pi) * numSlices_m);
    Px_m = Px / Physics::c;
    Py_m = Py / Physics::c;
    // dx02_m = dx0_m;
    // dy02_m = dy0_m;
    // dfi_x_m = 180.0 * dfi_x / Physics::pi;
    // dfi_y_m = 180.0 * dfi_y / Physics::pi;
    Ez_m = 1.0e-6 * Efz;
    Bz_m = Bfz;

    //in [mrad]
    emtn_m = Vector_t(ex, ey, 0.0);
    norm_emtn_m = Vector_t(nex, ney, 0.0);

    maxX_m = Vector_t(RxMax_m, RyMax_m, zMax);
    maxP_m = Vector_t(PxMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, PyMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, PzMax);

    //minX in the T-T sense is -RxMax_m
    minX_m = Vector_t(-RxMax_m, -RyMax_m, zMin);
    minP_m = Vector_t(-PxMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, -PyMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, PzMin);
    //minP_m = Vector_t(-maxP_m[0], -maxP_m[1], PzMin);

    /* PxRms is the rms of the divergence in x direction in [rad]: x'
     * x' = Px/P0 -> Px = x' * P0 = x' * \beta/c*E (E is the particle total
     * energy)
     * Pz* in [\beta\gamma]
     * \pi from [rad]
     */
    sigmax_m = Vector_t(RxRms / 2.0, RyRms / 2.0, zRms);
    sigmap_m = Vector_t(PxRms * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi / 2.0, PyRms * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi / 2.0, PzRms);

    dEdt_m = PzRms * mc2e * b0;

    IpplTimings::stopTimer(statParamTimer_m);
}

void EnvelopeBunch::runStats(EnvelopeBunchParameter sp, double *xAvg, double *xMax, double *xMin, double *rms, int *nValid) {
    int i, nV = 0;
    std::vector<double> v(numMySlices_m, 0.0);

    //FIXME: why from 1 to n-1??
    switch(sp) {
        case sp_beta:      // normalized velocity (total) [-]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_beta];
            }
            break;
        case sp_gamma:     // Lorenz factor
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->computeGamma();
            }
            break;
        case sp_z:         // slice position [m]
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_z];
            }
            break;
        case sp_I:         // slice position [m]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_beta] > BETA_MIN1) && ((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()))
                    v[nV++] = 0.0; //FIXME: currentProfile_m->get(slices_m[i]->p[SLI_z], itype_lin);
            }
            break;
        case sp_Rx:        // beam size x [m]
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = 2.* slices_m[i]->p[SLI_x];
            }
            break;
        case sp_Ry:        // beam size y [m]
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = 2. * slices_m[i]->p[SLI_y];
            }
            break;
        case sp_Px:        // beam divergence x
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_px];
            }
            break;
        case sp_Py:        // beam divergence y
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_py];
            }
            break;
        case sp_Pz:
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_beta] * slices_m[i]->computeGamma();
            }
            break;
        case sp_x0:        // position centroid x [m]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_x0];
            }
            break;
        case sp_y0:        // position centroid y [m]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_y0];
            }
            break;
        case sp_px0:       // angular deflection centroid x
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_px0];
            }
            break;
        case sp_py0:      // angular deflection centroid y
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) v[nV++] = slices_m[i]->p[SLI_py0];
            }
            break;
        default :
            throw OpalException("EnvelopeBunch", "EnvelopeBunch::runStats() Undefined label (programming error)");
            break;
    }

    int nVTot = nV;
    MPI_Allreduce(MPI_IN_PLACE, &nVTot, 1, MPI_INT, MPI_SUM, Ippl::getComm());
    if(nVTot <= 0) {
        *xAvg = 0.0;
        *xMax = 0.0;
        *xMin = 0.0;
        *rms = 0.0;
        *nValid = 0;
    } else {
        double M1 = 0.0;
        double M2 = 0.0;
        double maxv = std::numeric_limits<double>::min();
        double minv = std::numeric_limits<double>::max();
        if(nV > 0) {
            M1 = v[0];
            M2 = v[0] * v[0];
            maxv = v[0];
            minv = v[0];
            for(i = 1; i < nV; i++) {
                M1 += v[i];
                M2 += v[i] * v[i];
                maxv = std::max(maxv, v[i]);
                minv = std::min(minv, v[i]);
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &M1, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
        MPI_Allreduce(MPI_IN_PLACE, &M2, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
        MPI_Allreduce(MPI_IN_PLACE, &maxv, 1, MPI_DOUBLE, MPI_MAX, Ippl::getComm());
        MPI_Allreduce(MPI_IN_PLACE, &minv, 1, MPI_DOUBLE, MPI_MIN, Ippl::getComm());

        *xAvg = M1 / nVTot;
        *xMax = maxv;
        *xMin = minv;

        //XXX: ok so this causes problems. E.g in the case of all transversal
        //components we cannot compare a rms_rx with rms_x (particles)
        //produced by the T-Tracker

        //in case of transversal stats we want to calculate the rms
        //*rms = sqrt(M2/nVTot);

        //else a sigma
        //*sigma = sqrt(M2/nVTot - M1*M1/(nVTot*nVTot));
        //this is sigma:
        //*rms = sqrt(M2/nVTot - M1*M1/(nVTot*nVTot));
        *nValid = nVTot;

        if(sp == sp_Rx ||
           sp == sp_Ry ||
           sp == sp_Px ||
           sp == sp_Py)
            *rms = sqrt(M2 / nVTot); //2.0: half of the particles
        else
            *rms = sqrt(M2 / nVTot - M1 * M1 / (nVTot * nVTot));
    }
}

void EnvelopeBunch::calcEmittance(double *emtnx, double *emtny, double *emtx, double *emty, int *nValid) {
    double sx = 0.0;
    double sxp = 0.0;
    double sxxp = 0.0;
    double sy = 0.0;
    double syp = 0.0;
    double syyp = 0.0;
    double betagamma = 0.0;

    // find the amount of active slices
    int nV = 0;
    for(int i = 0; i < numMySlices_m; i++) {
        if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid())
            nV++;
    }

    if(nV > 0) {
        int i1 = nV;
        nV = 0;
        double bg = 0.0;
        for(int i = 0; i < i1; i++) {
            if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) {
                nV++;

                if(solver_m & sv_radial) {
                    assert(i < numMySlices_m);
                    auto sp = slices_m[i];
                    bg = sp->p[SLI_beta] * sp->computeGamma();

                    double pbc = bg * sp->p[SLI_px] / (sp->p[SLI_beta] * Physics::c);
                    sx   += sp->p[SLI_x] * sp->p[SLI_x];
                    sxp  += pbc * pbc;
                    sxxp += sp->p[SLI_x] * pbc;

                    pbc   = bg * sp->p[SLI_py] / (sp->p[SLI_beta] * Physics::c);
                    sy   += sp->p[SLI_y] * sp->p[SLI_y];
                    syp  += pbc * pbc;
                    syyp += sp->p[SLI_y] * pbc;

                    betagamma += sqrt(1 + sp->p[SLI_px] * sp->p[SLI_px] + sp->p[SLI_py] * sp->p[SLI_py]);
                }
            }
        }
    }

    int nVToT = nV;
    reduce(nVToT, nVToT, OpAddAssign());
    if(nVToT == 0) {
        *emtnx = 0.0;
        *emtny = 0.0;
        *emtx = 0.0;
        *emty = 0.0;
        *nValid = 0;
    } else {
        reduce(sx, sx, OpAddAssign());
        reduce(sy, sy, OpAddAssign());
        reduce(sxp, sxp, OpAddAssign());
        reduce(syp, syp, OpAddAssign());
        reduce(sxxp, sxxp, OpAddAssign());
        reduce(syyp, syyp, OpAddAssign());
        reduce(betagamma, betagamma, OpAddAssign());
        sx /= nVToT;
        sy /= nVToT;
        sxp /= nVToT;
        syp /= nVToT;
        sxxp /= nVToT;
        syyp /= nVToT;

        *emtnx = sqrt(sx * sxp - sxxp * sxxp + emtnx0_m * emtnx0_m + emtbx0_m * emtbx0_m);
        *emtny = sqrt(sy * syp - syyp * syyp + emtny0_m * emtny0_m + emtby0_m * emtby0_m);

        betagamma /= nVToT;
        betagamma *= sqrt(1.0 - (1 / betagamma) * (1 / betagamma));
        *emtx = *emtnx / betagamma;
        *emty = *emtny / betagamma;
        *nValid = nVToT;
    }
}

void EnvelopeBunch::calcEnergyChirp(double *g0, double *dgdt, double *gInc, int *nValid) {
    std::vector<double> dtl(numMySlices_m, 0.0);
    std::vector<double> b(numMySlices_m, 0.0);
    std::vector<double> g(numMySlices_m, 0.0);

    // defaults
    *g0   = 1.0;
    *dgdt = 0.0;
    *gInc = 0.0;

    double zAvg = 0.0;
    double gAvg = 0.0;
    int j = 0;
    for(int i = 0; i < numMySlices_m; i++) {
        std::shared_ptr<EnvelopeSlice> cs = slices_m[i];
        if((cs->is_valid()) && (cs->p[SLI_z] > zCat_m)) {
            zAvg    += cs->p[SLI_z];
            dtl[i-j]  = cs->p[SLI_z];
            b[i-j]   = cs->p[SLI_beta];
            g[i-j]   = cs->computeGamma();
            gAvg    += g[i-j];
        } else
            ++j;
    }
    int nV = numMySlices_m - j;

    int nVTot = nV;
    reduce(nVTot, nVTot, OpAddAssign());
    if(nVTot > 0) {
        reduce(gAvg, gAvg, OpAddAssign());
        reduce(zAvg, zAvg, OpAddAssign());
        gAvg = gAvg / nVTot;
        zAvg = zAvg / nVTot;
        *g0  = gAvg;
    }

    // FIXME: working with global arrays
    // make dtl, b and g global
    if(nVTot > 2) {
        std::vector<double> dtG(nVTot, 0.0);
        std::vector<double> bG(nVTot, 0.0);
        std::vector<double> gG(nVTot, 0.0);

        int numproc = Ippl::Comm->getNodes();
        std::vector<int> numsend(numproc, nV);
        std::vector<int> offsets(numproc);
        std::vector<int> offsetsG(numproc);
        std::vector<int> zeros(numproc, 0);

        MPI_Allgather(&nV, 1, MPI_INT, &offsets.front(), 1, MPI_INT, Ippl::getComm());
        offsetsG[0] = 0;
        for(int i = 1; i < numproc; i++) {
            if(offsets[i-1] == 0)
                offsetsG[i] = 0;
            else
                offsetsG[i] = offsetsG[i-1] + offsets[i-1];
        }

        MPI_Alltoallv(&dtl.front(), &numsend.front(), &zeros.front(), MPI_DOUBLE,
                      &dtG.front(), &offsets.front(), &offsetsG.front(), MPI_DOUBLE,
                      Ippl::getComm());

        MPI_Alltoallv(&b.front(), &numsend.front(), &zeros.front(), MPI_DOUBLE,
                      &bG.front(), &offsets.front(), &offsetsG.front(), MPI_DOUBLE,
                      Ippl::getComm());

        MPI_Alltoallv(&g.front(), &numsend.front(), &zeros.front(), MPI_DOUBLE,
                      &gG.front(), &offsets.front(), &offsetsG.front(), MPI_DOUBLE,
                      Ippl::getComm());

        double dum2, dum3, dum4, rms, gZero, gt;

        // convert z to t
        for(int i = 0; i < nVTot; i++) {
            dtG[i] = (dtG[i] - zAvg) / (bG[i] * Physics::c);
        }

        // chrip and uncorrelated energy sread
        linfit(&dtG[0], &gG[0], nVTot, &gZero, &gt, &dum2, &dum3, &dum4);
        *dgdt = gt;

        rms = 0.0;
        for(int i = 0; i < nVTot; i++) {
            rms += pow(gG[i] - gZero - gt * dtG[i], 2);
        }

        *gInc = sqrt(rms / nVTot);
    }
}


void EnvelopeBunch::distributeSlices(int nSlice) {
    numSlices_m = nSlice;
    int rank = Ippl::Comm->myNode();
    int numproc = Ippl::Comm->getNodes();
    numMySlices_m = nSlice / numproc;
    if(numMySlices_m < 13) {
        if(rank == 0) {
            numMySlices_m = 14;
        } else {
            numMySlices_m = (nSlice - 14) / (numproc - 1);
            if(rank - 1 < (nSlice - 14) % (numproc - 1))
                numMySlices_m++;
        }
    } else {
        if(rank < nSlice % numproc)
            numMySlices_m++;
    }

    mySliceStartOffset_m = rank * ((int)numSlices_m / numproc);
    if(rank < numSlices_m % numproc)
        mySliceStartOffset_m += rank;
    else
        mySliceStartOffset_m += numSlices_m % numproc;
    mySliceEndOffset_m = mySliceStartOffset_m + numMySlices_m - 1;
}


void EnvelopeBunch::createBunch() {

    // Memory issue when going to larger number of processors (this is
    // also the reason for the 30 slices hard limit per core!)
    // using heuristic (until problem is located)
    //size_t nSlices = (3 / 100.0 + 1.0) * numMySlices_m;

    size_t nSlices = getLocalNum();

    KR = std::unique_ptr<Vector_t[]>(new Vector_t[nSlices]);
    KT = std::unique_ptr<Vector_t[]>(new Vector_t[nSlices]);
    EF = std::unique_ptr<Vector_t[]>(new Vector_t[nSlices]);
    BF = std::unique_ptr<Vector_t[]>(new Vector_t[nSlices]);

    z_m.resize(numSlices_m, 0.0);
    b_m.resize(numSlices_m, 0.0);

    for(unsigned int i = 0; i < nSlices; i++) {
        KR[i] = Vector_t(0);
        KT[i] = Vector_t(0);
        BF[i] = Vector_t(0);
        EF[i] = Vector_t(0);
    }

    // set default DE solver method
    solver_m = sv_radial | sv_offaxis | sv_lwakes | sv_twakes;

    //FIXME: WHY?
    if(numSlices_m < 14)
        throw OpalException("EnvelopeBunch::createSlices", "use more than 13 slices");

    // defaults:
    Q_m        = 0.0;     // no charge
    t_m        = 0.0;     // t = 0 s
    t_offset_m = 0.0;     // offset time by tReset function
    emtnx0_m   = 0.0;
    emtny0_m   = 0.0;     // no intrinsic emittance
    emtbx0_m   = 0.0;
    emtby0_m   = 0.0;     // no intrinsic emittance Bush effect
    Bz0_m      = 0.0;     // no magnetic field on creation (cathode)
    dx0_m      = 0.0;
    dy0_m      = 0.0;     // no offset of coordinate system
    dfi_x_m    = 0.0;
    dfi_y_m    = 0.0;     // no rotation of coordinate system

    slices_m.resize(nSlices);
    for( auto & slice : slices_m ) {
        slice = std::shared_ptr<EnvelopeSlice>(new EnvelopeSlice());
    }
    //XXX: not supported by clang at the moment
    //std::generate(s.begin(), s.end(),
        //[]()
        //{
            //return std::shared_ptr<EnvelopeSlice>(new EnvelopeSlice());
        //});

    Esct_m.resize(nSlices, 0.0);
    G_m.resize(nSlices, 0.0);
    Exw_m.resize(nSlices, 0.0);
    Eyw_m.resize(nSlices, 0.0);
    Ezw_m.resize(nSlices, 0.0);

    for(unsigned int i = 0; i < nSlices; i++) {
        Esct_m[i] = 0.0;
        G_m[i]    = 0.0;
        Ezw_m[i]  = 0.0;
        Exw_m[i]  = 0.0;
        Eyw_m[i]  = 0.0;
    }

    currentProfile_m = NULL;
    I0avg_m = 0.0;
    dStat_m = ds_fieldsSynchronized | ds_slicesSynchronized;
}

//XXX: convention: s[n] is the slice closest to the cathode and all positions
//are negative (slices left of cathode)
//XXX: every processor fills its slices in bins (every proc holds all bins,
//some empty)
void EnvelopeBunch::setBinnedLShape(EnvelopeBunchShape shape, double z0, double emission_time_s, double frac) {
    unsigned int n2 = numSlices_m / 2;
    double sqr2 = sqrt(2.0), v;
    double bunch_width = 0.0;

    switch(shape) {
        case bsRect:

            bunch_width = Physics::c * emission_time_s * slices_m[0]->p[SLI_beta];
//            std::cout << "bunch_width = " << bunch_width << " SLI_beta= " << slices_m[0]->p[SLI_beta] << std::endl;
            for(int i = 0; i < numMySlices_m; i++) {
                slices_m[i]->p[SLI_z] = -(((numSlices_m - 1) - (mySliceStartOffset_m + i)) * bunch_width) / numSlices_m;
            }
            I0avg_m = Q_m * Physics::c / fabs(2.0 * emission_time_s);
            break;

        case bsGauss:
            if(n2 >= mySliceStartOffset_m && n2 <= mySliceEndOffset_m)
                slices_m[n2-mySliceStartOffset_m]->p[SLI_z] = z0;

            for(int i = 1; i <= numSlices_m / 2; i++) {
                rootValue = 1.0 - 2.0 * i * frac / (numSlices_m + 1);
                v = fabs(emission_time_s) * sqr2 * findRoot(erfRoot, 0.0, 5.0, 1.0e-5) * (emission_time_s < 0.0 ? Physics::c : 1.0);

                if((n2 + i) >= mySliceStartOffset_m && (n2 + i) <= mySliceEndOffset_m)
                    slices_m[n2+i-mySliceStartOffset_m]->p[SLI_z] = z0 + v * slices_m[n2+i-mySliceStartOffset_m]->p[SLI_beta];

                if((n2 - i) >= mySliceStartOffset_m && (n2 - i) <= mySliceEndOffset_m)
                    slices_m[n2-i-mySliceStartOffset_m]->p[SLI_z] = z0 - v * slices_m[n2-i-mySliceStartOffset_m]->p[SLI_beta];
            }

            I0avg_m = 0.0;
            break;
    }

    double gz0 = 0.0, gzN = 0.0;
    if(Ippl::Comm->myNode() == 0)
        gz0 = slices_m[0]->p[SLI_z];
    if(Ippl::Comm->myNode() == Ippl::Comm->getNodes() - 1)
        gzN = slices_m[numMySlices_m-1]->p[SLI_z];
    reduce(gz0, gz0, OpAddAssign());
    reduce(gzN, gzN, OpAddAssign());

    hbin_m = (gzN - gz0) / nebin_m;

    // initialize all bins with an empty vector
    for(int i = 0; i < nebin_m; i++) {
        std::vector<int> tmp;
        this->bins_m.push_back(tmp);
    }

    // add all slices to appropriated bins
    int bin_i = 0, slice_i = 0;
    while(slice_i < numMySlices_m) {
        if((bin_i + 1) * hbin_m < slices_m[slice_i]->p[SLI_z] - gz0) {
            bin_i++;
        } else {
            this->bins_m[bin_i].push_back(slice_i);
            slice_i++;
        }
    }

    ////XXX: debug output
    /*
    for(unsigned int i = 0; i < bins_m.size(); i++) {
      if(bins_m[i].size() > 0) {
    std::cout << Ippl::Comm->myNode() << ": Bin " << i << ": ";
    for(unsigned int j = 0; j < bins_m[i].size(); j++)
      std::cout << " " << mySliceStartOffset_m + bins_m[i][j] << "(" << slices_m[j]->p[SLI_z] << ")";
    std::cout << std::endl;
      }
    }
    */
    //find first bin containing slices
    firstBinWithValue_m = 0;
    for(; firstBinWithValue_m < nebin_m; firstBinWithValue_m++)
        if(bins_m[firstBinWithValue_m].size() != 0) break;

    backup();
}

void EnvelopeBunch::setTShape(double enx, double eny, double rx, double ry, double b0) {
  /*
    Inform msg("setTshape");
    msg << "set SLI_x to " << rx / 2.0 << endl;
    msg << "set SLI_y to " << ry / 2.0 << endl;
  */
    // set emittances
    emtnx0_m = enx;
    emtny0_m = eny;
    //FIXME: is rx here correct?
    emtbx0_m = Physics::q_e * rx * rx * Bz0_m / (8.0 * Physics::EMASS * Physics::c);
    emtby0_m = Physics::q_e * ry * ry * Bz0_m / (8.0 * Physics::EMASS * Physics::c);
    //msg << "set emtbx0 to " << emtbx0 << endl;
    //msg << "set emtby0 to " << emtby0 << endl;

    //XXX: rx = 2*rms_x
    for(int i = 0; i < numMySlices_m; i++) {
        slices_m[i]->p[SLI_x] = rx / 2.0; //rms_x
        slices_m[i]->p[SLI_y] = ry / 2.0; //rms_y
        slices_m[i]->p[SLI_px] = 0.0;
        slices_m[i]->p[SLI_py] = 0.0;
    }

    backup();
}

void EnvelopeBunch::setTOffset(double x0, double px0, double y0, double py0) {
    for(int i = 0; i < numMySlices_m; i++) {
        slices_m[i]->p[SLI_x0] = x0;
        slices_m[i]->p[SLI_px0] = px0;
        slices_m[i]->p[SLI_y0] = y0;
        slices_m[i]->p[SLI_py0] = py0;
    }
}

void EnvelopeBunch::setEnergy(double E0, double dE) {
    double g0    = 1.0 + (Physics::q_e * E0 / (Physics::EMASS * Physics::c * Physics::c));
    double dg    = fabs(dE) * Physics::q_e / (Physics::EMASS * Physics::c * Physics::c);
    double z0    = zAvg();

    for(int i = 0; i < numMySlices_m; i++) {
        double g = g0 + (slices_m[i]->p[SLI_z] - z0) * dg;
        slices_m[i]->p[SLI_beta] = sqrt(1.0 - (1.0 / (g * g)));
    }

    backup();
}

void EnvelopeBunch::synchronizeSlices() {

    for(int i = 0; i < numSlices_m; i++) {
        z_m[i] = 0.0;
        b_m[i] = 0.0;
    }
    for(int i = 0; i < numMySlices_m; i++) {
        b_m[mySliceStartOffset_m+i] = slices_m[i]->p[SLI_beta];
        z_m[mySliceStartOffset_m+i] = slices_m[i]->p[SLI_z];
    }

    MPI_Allreduce(MPI_IN_PLACE, &(z_m[0]), numSlices_m, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    MPI_Allreduce(MPI_IN_PLACE, &(b_m[0]), numSlices_m, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
}

void EnvelopeBunch::calcI() {
    Inform msg("calcI ");

    static int already_called = 0;
    if((dStat_m & ds_currentCalculated) || (already_called && (Q_m <= 0.0)))
        return;
    already_called = 1;

    std::vector<double> z1(numSlices_m, 0.0);
    std::vector<double> b (numSlices_m, 0.0);
    double bSum = 0.0;
    double dz2Sum = 0.0;
    int n1 = 0;
    int my_start = 0, my_end = 0;

    for(int i = 0; i < numSlices_m; i++) {
        if(b_m[i] > 0.0) {
            if((unsigned int) i == mySliceStartOffset_m)
                my_start = n1;
            if((unsigned int) i == mySliceEndOffset_m)
                my_end = n1;

            b[n1] = b_m[i];
            z1[n1] = z_m[i];
            if(n1 > 0)
                dz2Sum += ((z1[n1] - z1[n1-1]) * (z1[n1] - z1[n1-1]));
            bSum += b_m[i];
            n1++;
        }
    }
    if(Ippl::Comm->myNode() == 0)
        my_start++;

    if(n1 < 2) {
        msg << "n1 (= " << n1 << ") < 2" << endl;
        //throw OpalException("EnvelopeBunch", "Insufficient points to calculate the current (n1)");
        currentProfile_m = std::unique_ptr<Profile>(new Profile(0.0));
        return;
    }

    double sigma_dz = sqrt(dz2Sum / (n1 - 1));
    double beta = bSum / n1;

    //sort beta's according to z1 with temporary vector of pairs
    std::vector<std::pair<double,double>> z1_b(numSlices_m);
    for (int i=0; i<numSlices_m; i++) {
        z1_b[i] = std::make_pair(z1[i],b[i]);
    }

    std::sort(z1_b.begin(), z1_b.end(),
              [](const std::pair<double,double> &left, const std::pair<double,double> &right)
              {return left.first < right.first;});

    for (int i=0; i<numSlices_m; i++) {
        z1[i] = z1_b[i].first;
        b [i] = z1_b[i].second;
    }

    double q = Q_m > 0.0 ? Q_m / numSlices_m : Physics::q_e;
    double dz = 0.0;

    // 1. Determine current from the slice distance
    std::vector<double> I1(n1, 0.0);

    //limit the max current to 5x the sigma value to reduce noise problems
    double dzMin = 0.2 * sigma_dz;

    //FIXME: Instead of using such a complicated mechanism (store same value
    //again) to enforce only to use a subsample of all points, use a vector to
    //which all values get pushed, no need for elimination later.

    int vend = my_end;
    if(Ippl::Comm->myNode() == Ippl::Comm->getNodes() - 1)
        vend--;

    //XXX: THIS LOOP IS THE EXPENSIVE PART!!!
    int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,dz)
#endif //_OPENMP
    for(i = my_start; i <= vend; i++) {
        j = 1;
        do {
            dz = fabs(z1[i+j] - z1[i-j]);
            j++;
        } while((dz < dzMin * (j - 1)) && ((i + j) < n1) && ((i - j) >= 0));

        if((dz >= dzMin * (j - 1)) && ((i + j) < n1) && ((i - j) >= 0)) {
            I1[i] = 0.25 * q * Physics::c * (b[i+j] + b[i-j]) / (dz * (j - 1));   // WHY 0.25?
        } else {
            I1[i] = 0.0;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &(I1[0]), n1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    for(int i = 1; i < n1 - 1; i++) {
        if(I1[i] == 0.0)
            I1[i] = I1[i-1];
    }
    I1[0]    = I1[1];
    I1[n1-1] = I1[n1-2];


    //2. Remove points with identical z-value and then smooth the current profile
    double zMin = zTail();
    double zMax = zHead();
    dz = (zMax - zMin) / numSlices_m; // create a window of the average slice distance
    std::vector<double> z2(n1, 0.0);
    std::vector<double> I2(n1, 0.0);
    double Mz1 = 0.0;
    double MI1 = 0.0;
    int np = 0;
    j = 0;

    // XXX: COMPUTE ON ALL NODES
    // first value
    while((j < n1) && ((z1[j] - z1[0]) <= dz)) {
        Mz1 += z1[j];
        MI1 += I1[j];
        ++j;
        ++np;
    }
    z2[0] = Mz1 / np;
    I2[0] = MI1 / np;

    // following values
    int k = 0;
    for(int i = 1; i < n1; i++) {
        //for(int i = my_start; i <= my_end; i++) {
        // add new points
        int j = 0;
        while(((i + j) < n1) && ((z1[i+j] - z1[i]) <= dz)) {
            if((z1[i+j] - z1[i-1]) > dz) {
                Mz1 += z1[i+j];
                MI1 += I1[i+j];
                ++np;
            }
            ++j;
        }

        // remove obsolete points
        j = 1;
        while(((i - j) >= 0) && ((z1[i-1] - z1[i-j]) < dz)) {
            if((z1[i] - z1[i-j]) > dz) {
                Mz1 -= z1[i-j];
                MI1 -= I1[i-j];
                --np;
            }
            ++j;
        }
        //z2_temp[i] = Mz1 / np;
        //I2_temp[i] = MI1 / np;
        z2[i-k] = Mz1 / np;
        I2[i-k] = MI1 / np;

        // make sure there are no duplicate z coordinates
        if(z2[i-k] <= z2[i-k-1]) {
            I2[i-k-1] = 0.5 * (I2[i-k] + I2[i-k-1]);
            k++;
        }
    }

    //MPI_Allreduce(MPI_IN_PLACE, &(z2_temp[0]), n1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    //MPI_Allreduce(MPI_IN_PLACE, &(I2_temp[0]), n1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());

    ////FIXME: we dont need copy of z2 and I2: z2[i-k] = z2[i];
    //int k = 0;
    //for(int i = 1; i < n1; i++) {
    //z2[i-k] = z2_temp[i];
    //I2[i-k] = I2_temp[i];
    //// make sure there are no duplicate z coordinates
    //if(z2[i-k] <= z2[i-k-1]) {
    //I2[i-k-1] = 0.5 * (I2[i-k] + I2[i-k-1]);
    //k++;
    //}
    //}

    //free(z2_temp);
    //free(I2_temp);

    int n2 = n1 - k;
    if(n2 < 1) {
        msg << "Insufficient points to calculate the current (m = " << n2 << ")" << endl;
        currentProfile_m = std::unique_ptr<Profile>(new Profile(0.0));
    } else {
        // 3. smooth further
        if(n2 > 40) {
            sgSmooth(&I2.front(), n2, n2 / 20, n2 / 20, 0, 1);
        }

        // 4. create current profile
        currentProfile_m = std::unique_ptr<Profile>(new Profile(&z2.front(),
                                                    &I2.front(), n2));

        /**5. Normalize profile to match bunch charge as a constant
         * However, only normalize for sufficient beam energy
         */
        thisBunch = this;

        double Qcalc = 0.0;
        double z = zMin;
        dz = (zMax - zMin) / 99;
        for(int i = 1; i < 100; i++) {
            Qcalc += currentProfile_m->get(z, itype_lin);
            z += dz;
        }
        Qcalc *= (dz / (beta * Physics::c));
        currentProfile_m->scale((Q_m > 0.0 ? Q_m : Physics::q_e) / Qcalc);
    }

    dStat_m |= ds_currentCalculated;
}

void EnvelopeBunch::cSpaceCharge() {
    Inform msg("cSpaceCharge");

#ifdef USE_HOMDYN_SC_MODEL
    int icON = 1;
    double zhead = zHead();
    double ztail = zTail();
    double L = zhead - ztail;

    for(int i = 0; i < numMySlices_m; i++) {

        if(slices_m[i]->p[SLI_z] > 0.0)  {
            double zeta = slices_m[i]->p[SLI_z] - ztail;
            double xi   = slices_m[i]->p[SLI_z] + zhead;
            double sigma_av = (slices_m[i]->p[SLI_x] + slices_m[i]->p[SLI_y]) / 2;
            double R = 2 * sigma_av;
            double A = R / L / getGamma(i);

            double H1 = sqrt((1 - zeta / L) * (1 - zeta / L) + A * A) - sqrt((zeta / L) * (zeta / L) + A * A) - fabs(1 - zeta / L) + fabs(zeta / L);
            double H2 = sqrt((1 - xi / L) * (1 - xi / L) + A * A) - sqrt((xi / L) * (xi / L) + A * A) - fabs(1 - xi / L) + fabs(xi / L);

            //FIXME: Q_act or Q?
            //double Q = activeSlices_m * Q_m / numSlices_m;
            Esct_m[i] = (Q_m / 2 / Physics::pi / Physics::epsilon_0 / R / R) * (H1 - icON * H2);
            double G1 = (1 - zeta / L) / sqrt((1 - zeta / L) * (1 - zeta / L) + A * A) + (zeta / L) / sqrt((zeta / L) * (zeta / L) + A * A);
            double G2 = (1 - xi / L) / sqrt((1 - xi / L) * (1 - xi / L) + A * A) + (xi / L) / sqrt((xi / L) * (xi / L) + A * A);

            G_m[i] = (1 - getBeta(i) * getBeta(i)) * G1 - icON * (1 + getBeta(i) * getBeta(i)) * G2;
        }
    }
#else
    if(numSlices_m < 2) {
        msg << "called with insufficient slices (" << numSlices_m << ")" << endl;
        return;
    }

    if((Q_m <= 0.0) || (currentProfile_m->max() <= 0.0)) {
        msg << "Q or I_max <= 0.0, aborting calculation" << endl;
        return;
    }

    for(int i = 0; i < numMySlices_m; ++i) {
        Esct[i] = 0.0;
        G[i]    = 0.0;
    }

    double Imax = currentProfile_m->max();
    int nV = 0;
    double sm = 0.0;
    double A0 = 0.0;
    std::vector<double> xi(numSlices_m, 0.0);

    for(int i = 0; i < numMySlices_m; i++) {
        if(slices_m[i]->p[SLI_z] > 0.0) {
            nV++;
            A0 = 4.0 * slices_m[i]->p[SLI_x] * slices_m[i]->p[SLI_y];
            sm += A0;
            xi[i+mySliceStartOffset_m] = A0 * (1.0 - slices_m[i]->p[SLI_beta] * slices_m[i]->p[SLI_beta]); // g2
        }
    }

    int nVTot = nV;
    MPI_Allreduce(MPI_IN_PLACE, &nVTot, 1, MPI_INT, MPI_SUM, Ippl::getComm());
    if(nVTot < 2) {
        msg << "Exiting, to few nV slices" << endl;
        return;
    }

    MPI_Allreduce(MPI_IN_PLACE, &xi[0], numSlices_m, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    MPI_Allreduce(MPI_IN_PLACE, &sm, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    A0 = sm / nVTot;
    double dzMin = 5.0 * Physics::c * Q_m / (Imax * numSlices_m);

    for(int localIdx = 0; localIdx < numMySlices_m; localIdx++) {
        double z0 = z_m[localIdx + mySliceStartOffset_m];
        if(z0 > 0.0 && slices_m[localIdx]->p[SLI_beta] > BETA_MIN1) {

            sm = 0.0;

            for(int j = 0; j < numSlices_m; j++) {
                double zj = z_m[j];
                double dz = fabs(zj - z0);
                if((dz > dzMin) && (zj > zCat)) {
                    double Aj = xi[j] / (dz * dz);
                    double v  = 1.0 - (1.0 / sqrt(1.0 + Aj));
                    if(zj > z0) {
                        sm -= v;
                    } else {
                        sm += v;
                    }
                }
            }

            // longitudinal effect
            double bt = slices_m[localIdx]->p[SLI_beta];
            double btsq = (bt - BETA_MIN1) / (BETA_MIN2 - BETA_MIN1) * (bt - BETA_MIN1) / (BETA_MIN2 - BETA_MIN1);
            Esct[localIdx] = (bt < BETA_MIN1 ? 0.0 : bt < BETA_MIN2 ? btsq : 1.0);
            Esct[localIdx] *= Q_m * sm / (Physics::two_pi * Physics::epsilon_0 * A0 * (nVTot - 1));
            G[localIdx] = currentProfile_m->get(z0, itype_lin) / Imax;

            // tweak to compensate for non-relativity
            if(bt < BETA_MIN2) {
                if(slices_m[localIdx]->p[SLI_beta] < BETA_MIN1)
                    G[localIdx] = 0.0;
                else
                    G[localIdx] *= btsq;
            }
        }
    }

    return;
#endif
}

double EnvelopeBunch::moveZ0(double zC) {
    zCat_m = zC;
    double dz = zC - zHead();
    if(dz > 0.0) {
        for(int i = 0; i < numMySlices_m; i++) {
            slices_m[i]->p[SLI_z] += dz;
        }
        backup(); // save the new state
        *gmsg << "EnvelopeBunch::moveZ0(): bunch moved with " << dz << " m to " << zCat_m << " m" << endl;
    }

    return dz;
}

double EnvelopeBunch::tReset(double dt) {
    double new_dt = dt;

    if(dt == 0.0) {
        new_dt = t_m;
        *gmsg << "EnvelopeBunch time reset at z = " << zAvg() << " m with: " << t_m << " s, new offset: " << t_m + t_offset_m << " s";
    }

    t_offset_m += new_dt;
    t_m        -= new_dt;

    return new_dt;
}

/** derivs for RK routine
 * Definition of equation set:
 *  ==========================
 Y[SLI_z]    = z      dz/dt  = beta*c*cos(a)
                      cos(a) = 1 - (px0^2 + py0^2)/c2
 Y[SLI_beta] = beta   db/dt  = (e0/mc)*(E_acc + E_sc)/gamma^3
 Y[SLI_x]    = x      dx/dt  = px = Y[SLI_px]
 Y[SLI_px]   = px     dpx/dt = f(x,beta) - (beta*gamma^2(db/dt)*px + Kr*x)
 Y[SLI_y]    = y      dy/dt  = py = Y[SLI_py]
 Y[SLI_py]   = py     dpy/dt = f(y,beta) - (beta*gamma^2(db/dt)*py + Kr*y)
 Y[SLI_x0]   = x0     dx0/dt = px0 = Y[SLI_px0]
 Y[SLI_px0]  = px0    dpx0/dt= -(beta*gamma^2(db/dt)*px0) + Kt*x0
 Y[SLI_y0]   = y0     dy0/dt = py0 = Y[SLI_py0]
 Y[SLI_py0]  = py0    dpy0/dt= -(beta*gamma^2(db/dt)*py0) + Kt*y0

 Transversal space charge blowup:
 f(x,beta) = c^2*I/(2Ia)/(x*beta*gamma^3)

ALT: SLI_z (commented by Rene)
    // dYdt[SLI_z] = Y[SLI_beta]*c*sqrt(1.0 - (pow(Y[SLI_px0],2) + pow(Y[SLI_py0],2))/pow(c*Y[SLI_beta],2));
    // dYdt[SLI_z] = Y[SLI_beta]*c*cos(Y[SLI_px0]/Y[SLI_beta]/c)*cos(Y[SLI_py0]/Y[SLI_beta]/c);

 */
void EnvelopeBunch::derivs(double tc, double Y[], double dYdt[]) {
    double g2 = 1.0 / (1.0 - Y[SLI_beta] * Y[SLI_beta]);
    double g  = sqrt(g2);
    double g3 = g2 * g;

    double alpha = sqrt(Y[SLI_px0] * Y[SLI_px0] + Y[SLI_py0] * Y[SLI_py0]) / Y[SLI_beta] / Physics::c;
    /// \f[ \dot{z} = \beta c cos(\alpha) \f]
    dYdt[SLI_z]  = Y[SLI_beta] * Physics::c * cos(alpha);

    //NOTE: (reason for - when using Esl(2)) In OPAL we use e_0 with a sign.
    // The same issue with K factors (accounted in K factor calculation)
    //
    /// \f[ \dot{\beta} = \frac{e_0}{mc\gamma^3} \left(E_{z,\text{ext}} + E_{z,\text{sc}} + E_{z, \text{wake}} \right)
    dYdt[SLI_beta] = Physics::e0mc * (-Esl_m(2) + Esct_m[currentSlice_m] + Ezw_m[currentSlice_m]) / g3;

    /// \f[ \beta * \gamma^2 * \dot{\beta} \f]
    double bg2dbdt = Y[SLI_beta] * g2 * dYdt[SLI_beta];

    if(solver_m & sv_radial) {
        /// minimum spot-size due to emittance
        /// \f[ \left(\frac{\epsilon_n c}{\gamma}\right)^2 \f]
        double enxc2 = pow((emtnx0_m + emtbx0_m) * Physics::c / (Y[SLI_beta] * g), 2);
        double enyc2 = pow((emtny0_m + emtby0_m) * Physics::c / (Y[SLI_beta] * g), 2);


#ifdef USE_HOMDYN_SC_MODEL
        //FIXME: Q_act or Q?
        //double Q = activeSlices_m * Q_m / numSlices_m;
        double kpc  = 0.5 * Physics::c * Physics::c * (Y[SLI_beta] * Physics::c) * activeSlices_m * Q_m / numSlices_m / (curZHead_m - curZTail_m) / Physics::Ia;

        /// \f[ \dot{\sigma} = p \f]
        dYdt[SLI_x] = Y[SLI_px];
        dYdt[SLI_y] = Y[SLI_py];

        double sigma_av = (Y[SLI_x] + Y[SLI_y]) / 2;

        /// \f[ \ddot{\sigma} = -\gamma^2\beta\dot{\beta}\dot{\sigma} - K\sigma + 2c^2\left(\frac{I}{2I_0}\right)\frac{G}{\beta R^2}(1-\beta)^2 \sigma + \left(\frac{\epsilon_n c}{\gamma}\right)^2 \frac{1}{\sigma^3} \f]
        dYdt[SLI_px] = -bg2dbdt * Y[SLI_px] - KRsl_m(0) * Y[SLI_x] + (kpc / sigma_av / Y[SLI_beta] / g / 2) * G_m[currentSlice_m] + enxc2 / g3;
        dYdt[SLI_py] = -bg2dbdt * Y[SLI_py] - KRsl_m(1) * Y[SLI_y] + (kpc / sigma_av / Y[SLI_beta] / g / 2) * G_m[currentSlice_m] + enyc2 / g3;
#else
        // transverse space charge
        // somewhat strange: I expected: c*c*I/(2*Ia) (R. Bakker)
        //XXX: Renes version, I[cs] already in G
        double kpc  = 0.5 * Physics::c * Physics::c * currentProfile_m->max() / Physics::Ia;

        /// \f[ \dot{\sigma} = p \f]
        dYdt[SLI_x]  = Y[SLI_px];
        dYdt[SLI_y]  = Y[SLI_py];
        dYdt[SLI_px] = (kpc * G_m[currentSlice_m] / (Y[SLI_x] * Y[SLI_beta] * g3)) + enxc2 / g3 - (KRsl_m(0) * Y[SLI_x]) - (bg2dbdt * Y[SLI_px]);
        dYdt[SLI_py] = (kpc * G_m[currentSlice_m] / (Y[SLI_y] * Y[SLI_beta] * g3)) + enyc2 / g3 - (KRsl_m(1) * Y[SLI_y]) - (bg2dbdt * Y[SLI_py]);
#endif
    } else {
        /// \f[ \dot{\sigma} = p \f]
        dYdt[SLI_x]  = Y[SLI_px];
        dYdt[SLI_y]  = Y[SLI_py];
        dYdt[SLI_px] = 0.0;
        dYdt[SLI_py] = 0.0;
    }

    if(solver_m & sv_offaxis) {
        dYdt[SLI_x0]  = Y[SLI_px0];
        dYdt[SLI_y0]  = Y[SLI_py0];
        dYdt[SLI_px0] = -KTsl_m(0) - (bg2dbdt * Y[SLI_px0]) + Physics::e0m * (g * Exw_m[currentSlice_m]);
        dYdt[SLI_py0] = -KTsl_m(1) - (bg2dbdt * Y[SLI_py0]) + Physics::e0m * (g * Eyw_m[currentSlice_m]);
    } else {
        dYdt[SLI_x0]  = Y[SLI_px0];
        dYdt[SLI_y0]  = Y[SLI_py0];
        dYdt[SLI_px0] = 0.0;
        dYdt[SLI_py0] = 0.0;
    }
}


void EnvelopeBunch::computeSpaceCharge() {
    IpplTimings::startTimer(selfFieldTimer_m);

    // Calculate the current profile
    if(Q_m > 0.0) {

        IpplTimings::startTimer(calcITimer_m);
#ifndef USE_HOMDYN_SC_MODEL
        //XXX: only used in BET-SC-MODEL
        // make sure all processors have all global betas and positions
        synchronizeSlices();
        calcI();
#endif
        IpplTimings::stopTimer(calcITimer_m);

        // the following assumes space-charges do not change significantly over nSc steps
        IpplTimings::startTimer(spaceChargeTimer_m);
        cSpaceCharge();
        IpplTimings::stopTimer(spaceChargeTimer_m);

        //if(numSlices_m > 40) {
        // smooth Esct to get rid of numerical noise
        //double Esc = 0;
        //sgSmooth(Esct,n,n/20,n/20,0,4);
        //}
    } else {
        currentProfile_m = std::unique_ptr<Profile>(new Profile(0.0));
    }

    IpplTimings::stopTimer(selfFieldTimer_m);
}


void EnvelopeBunch::timeStep(double tStep, double _zCat) {
    Inform msg("tStep");
    static int msgParsed = 0;

    // default accuracy of integration
    double eps = 1.0e-4;
    double time_step_s = tStep;

    thisBunch = this;
    zCat_m     = _zCat;

    // backup last stage before new execution
    backup();

    //FIXME: HAVE A PROBLEM WITH EMISSION (EMISSION OFF GIVES GOOD RESULTS)
    activeSlices_m = numSlices_m;

    //if(lastEmittedBin_m < nebin_m) {

    //time_step_s = emission_time_step_;
    //size_t nextBin = nebin_m - 1 - lastEmittedBin_m;
    //if(Ippl::Comm->myNode() == 0) cout << "Emitting bin " << nextBin << endl;

    //double gz0 = 0.0;
    //if(Ippl::Comm->myNode() == 0)
    //gz0 = slices_m[0]->p[SLI_z];
    //MPI_Allreduce(MPI_IN_PLACE, &gz0, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());

    ////XXX: bin holds local slice number
    //for(int j = 0; j < bins_m[nextBin].size(); j++) {
    //slices_m[bins_m[nextBin][j]]->p[SLI_z] -= gz0;
    //slices_m[bins_m[nextBin][j]]->p[SLI_z] -= hbin_m * nextBin;
    //slices_m[bins_m[nextBin][j]]->emitted = true;
    //cout << "\tEmitting slice " << mySliceStartOffset_m + bins_m[nextBin][j] << " at " << slices_m[bins_m[nextBin][j]]->p[SLI_z] << endl;
    //}

    //lastEmittedBin_m++;
    //activeSlices_m += bins_m[nextBin].size();
    //}

    //activeSlices_m = min(activeSlices_m+1, numSlices_m);
    //if(activeSlices_m != numSlices_m)
    //time_step_s = emission_time_step_;
    //else
    //time_step_s = tStep;

    curZHead_m = zHead();
    curZTail_m = zTail();

    for(int i = 0; i < numMySlices_m; i++) {
        // make the current slice index global in this class
        currentSlice_m = i;
        std::shared_ptr<EnvelopeSlice> sp = slices_m[i];

        //if(i + mySliceStartOffset_m == numSlices_m - activeSlices_m) {
        //sp->emitted = true;
        //sp->p[SLI_z] = 0.0;
        //std::cout << "emitting " << i + mySliceStartOffset_m << std::endl;
        //}

        // Assign values of K for certain slices
        KRsl_m = KR[i];
        KTsl_m = KT[i];
        Esl_m = EF[i];
        Bsl_m = BF[i];

        // only for slices already emitted
        if(true /*sp->emitted*/) {
            // set default allowed error for integration
            double epsLocal = eps;
            // mark that the integration was not successful yet
            int ode_result = 1;

            while(ode_result == 1) {

                if(solver_m & sv_fixedStep) {
                    rk4(&(sp->p[SLI_z]), SLNPAR, t_m, time_step_s, Gderivs);
                    ode_result = 0;
                } else {
                    int nok, nbad;
                    activeSlice_m = i;
                    ode_result = odeint(&(sp->p[SLI_z]), SLNPAR, t_m, t_m + time_step_s, epsLocal, 0.1 * time_step_s, 0.0, &nok, &nbad, Gderivs);
                }

                if(ode_result != 0) {
                    // restore the backup
                    sp->restore();
                    epsLocal *= 10.0;
                }
            }

            if(ode_result == 1) {
                // use fixed step integration if dynamic fails
                rk4(&(sp->p[SLI_z]), SLNPAR, t_m, time_step_s, Gderivs);

                if(msgParsed < 2) {
                    msg << "EnvelopeBunch::run() Switched to fixed step RK rountine for solving of DE at slice " << i << endl;
                    msgParsed = 2;
                }
            } else if((epsLocal != eps) && (msgParsed == 0)) {
                msg << "EnvelopeBunch::run() integration accuracy relaxed to " << epsLocal << " for slice " << i << " (ONLY FIRST OCCURENCE MARKED!)" << endl;
                msgParsed = 1;
            }
        }

        if(slices_m[i]->check()) {
            msg << "Slice " << i << " no longer valid at z = " <<  slices_m[i]->p_old[SLI_z] << " beta = " << slices_m[i]->p_old[SLI_beta] << endl;
            msg << "Slice " << i << " no longer valid at z = " <<  slices_m[i]->p[SLI_z] << " beta = " << slices_m[i]->p[SLI_beta] << endl;
            isValid_m = false;
            return;
        }
    }
    // mark that slices might not be synchronized (and space charge accordingly)
    dStat_m &= (!(ds_slicesSynchronized | ds_spaceCharge));

    /// mark calling of this function + update vars
    t_m += time_step_s;

    /// subtract average orbit for when tracking along the s-axis
    if(solver_m & sv_s_path) {
        int nV = 0;
        double ga = 0.0, x0a = 0.0, px0a = 0.0, y0a = 0.0, py0a = 0.0;
        double beta, fi_x, fi_y;

        //FIXME: BET calculates only 80 %, OPAL doesn't ?

        for(int i = 0; i < numMySlices_m; i++) {
            std::shared_ptr<EnvelopeSlice> sp  = slices_m[i];
            if((sp->p[SLI_z] >= zCat_m) && sp->is_valid()) {
                ++nV;
                ga  += sp->computeGamma();
                x0a += sp->p[SLI_x0];
                y0a += sp->p[SLI_y0];
                px0a += sp->p[SLI_px0];
                py0a += sp->p[SLI_py0];
            }
        }

        int nVTot = nV;
        reduce(nVTot, nVTot, OpAddAssign());
        if(nVTot > 0) {
            if(nV > 0) {
                reduce(ga, ga, OpAddAssign());
                reduce(x0a, x0a, OpAddAssign());
                reduce(px0a, px0a, OpAddAssign());
                reduce(y0a, y0a, OpAddAssign());
                reduce(py0a, py0a, OpAddAssign());
            }
            ga  = ga / nVTot;
            x0a = x0a / nVTot;
            px0a = px0a / nVTot;
            y0a = y0a / nVTot;
            py0a = py0a / nVTot;
        } else {
            msg << "EnvelopeBunch::run() No valid slices to subtract average" << endl;
        }
        beta = sqrt(1.0 - (1.0 / pow(ga, 2)));
        fi_x = px0a / Physics::c / beta;
        fi_y = py0a / Physics::c / beta;

        dx0_m -= x0a;
        dy0_m -= y0a;
        dfi_x_m -= fi_x;
        dfi_y_m -= fi_y;
        for(int i = 0; i < numMySlices_m; i++) {
            std::shared_ptr<EnvelopeSlice> sp = slices_m[i];

            sp->p[SLI_x0] -= x0a;
            sp->p[SLI_y0] -= y0a;
            sp->p[SLI_px0] -= px0a;
            sp->p[SLI_py0] -= py0a;
            sp->p[SLI_z] += (sp->p[SLI_x0] * sin(fi_x) + sp->p[SLI_y0] * sin(fi_y));
        }
    }
}

void EnvelopeBunch::initialize(int num_slices, double Q, double energy, double width, double emission_time, double frac, double current, double bunch_center, double bX, double bY, double mX, double mY, double Bz0, int nbin) {

#ifdef USE_HOMDYN_SC_MODEL
    *gmsg << "* Using HOMDYN space-charge model" << endl;
#else
    *gmsg << "* Using BET space-charge model" << endl;
#endif
    distributeSlices(num_slices);
    createBunch();

    this->setCharge(Q);
    this->setEnergy(energy);

    //FIXME: how do I have to set this (only used in Gauss)
    bunch_center = -1 * emission_time / 2.0;

    //FIXME: pass values
    nebin_m = nbin;

    lastEmittedBin_m = 0;
    activeSlices_m = 0;
    emission_time_step_m = emission_time / nebin_m;

    this->setBinnedLShape(bsRect, bunch_center, emission_time, frac);
    this->setTShape(mX, mY, bX, bY, Bz0);

    //FIXME:
    // 12: on-axis, radial, default track all
    this->setSolverParameter(12);
}

double EnvelopeBunch::AvBField() {
    double bf = 0.0;
    for(int slice = 0; slice < numMySlices_m; slice++) {
        for(int i = 0; i < 3; i++) {
            bf += BF[slice](i);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &bf, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    return bf / numSlices_m;
}

double EnvelopeBunch::AvEField() {
    double ef = 0.0;
    for(int slice = 0; slice < numMySlices_m; slice++) {
        for(int i = 0; i < 3; i++) {
            ef += EF[slice](i);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &ef, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    return ef / numSlices_m;
}

double EnvelopeBunch::Eavg() {
    int nValid = 0;
    double sum = 0.0;
    for(int i = 0; i < numMySlices_m; i++) {
        if((slices_m[i]->p[SLI_z] > zCat_m) && slices_m[i]->is_valid()) {
            sum += slices_m[i]->computeGamma();
            nValid++;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &nValid, 1, MPI_INT, MPI_SUM, Ippl::getComm());
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    sum /= nValid;
    return (nValid > 0 ? ((Physics::EMASS * Physics::c * Physics::c / Physics::q_e) * (sum - 1.0)) : 0.0);
}

double EnvelopeBunch::get_sPos() {
    //FIXME: find reference position = centroid?
    double refpos = 0.0;
    size_t count = 0;
    for(int i = 0; i < numMySlices_m; i++) {
        if(slices_m[i]->p[SLI_z] > 0.0) {
            refpos += slices_m[i]->p[SLI_z];
            count++;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_LONG, MPI_SUM, Ippl::getComm());
    MPI_Allreduce(MPI_IN_PLACE, &refpos, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    return refpos / count;
}

double EnvelopeBunch::zAvg() {
    int nV = 0;
    double sum = 0.0;
    for(int i = 0; i < numMySlices_m; i++) {
        if(slices_m[i]->is_valid()) {
            sum += slices_m[i]->p[SLI_z];
            nV++;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &nV, 1, MPI_INT, MPI_SUM, Ippl::getComm());
    if(nV < 1) {
        isValid_m = false;
        return -1;
        //throw OpalException("EnvelopeBunch", "EnvelopeBunch::zAvg() no valid slices left");
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, Ippl::getComm());
    return (sum / nV);
}

double EnvelopeBunch::zTail() {
    double min;

    int i = 0;
    while((i < numMySlices_m) && (!slices_m[i]->is_valid()))
        i++;
    //throw OpalException("EnvelopeBunch", "EnvelopeBunch::zTail() no valid slices left");
    if(i == numMySlices_m)
        isValid_m = false;
    else
        min = slices_m[i]->p[SLI_z];

    for(i = i + 1; i < numMySlices_m; i++)
        if((slices_m[i]->p[SLI_z] < min) && (slices_m[i]->is_valid()))
            min = slices_m[i]->p[SLI_z];

    //reduce(min, min, OpMinAssign());
    MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, Ippl::getComm());
    return min;
}

double EnvelopeBunch::zHead() {
    double max;

    int i = 0;
    while((i < numMySlices_m) && (slices_m[i]->is_valid() == 0))
        i++;
    //throw OpalException("EnvelopeBunch", "EnvelopeBunch::zHead() no valid slices left");
    if(i == numMySlices_m)
        isValid_m = false;
    else
        max = slices_m[i]->p[SLI_z];

    for(i = 1; i < numMySlices_m; i++)
        if(slices_m[i]->p[SLI_z] > max) max = slices_m[i]->p[SLI_z];

    //reduce(max, max, OpMaxAssign());
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, Ippl::getComm());
    return max;
}


Inform &EnvelopeBunch::slprint(Inform &os) {
    if(this->getTotalNum() != 0) {  // to suppress Nan's
        os << "* ************** S L B U N C H ***************************************************** " << endl;
        os << "* NSlices= " << this->getTotalNum() << " Qtot= " << Q_m << endl; //" [nC]  Qi= " << std::abs(qi_m) << " [C]" << endl;
        os << "* Emean= " << get_meanKineticEnergy() * 1e-6 << " [MeV]" << endl;
        os << "* dT= " << this->getdT() << " [s]" << endl;
        os << "* spos= " << this->zAvg() << " [m]" << endl;

        //os << "* rmax= " << rmax_m << " [m]" << endl;
        //os << "* rmin= " << rmin_m << " [m]" << endl;
        //os << "* rms beam size= " << rrms_m << " [m]" << endl;
        //os << "* rms momenta= " << prms_m << " [beta gamma]" << endl;
        //os << "* mean position= " << rmean_m << " [m]" << endl;
        //os << "* mean momenta= " << pmean_m << " [beta gamma]" << endl;
        //os << "* rms emmitance= " << eps_m << " (not normalized)" << endl;
        //os << "* rms correlation= " << rprms_m << endl;
        //os << "* tEmission= " << getTEmission() << " [s]" << endl;
        os << "* ********************************************************************************** " << endl;

        //Inform::FmtFlags_t ff = os.flags();
        //os << scientific;
        //os << "* ************** B U N C H ********************************************************* " << endl;
        //os << "* NSlices= " << this->getTotalNum() << " Qtot= " << abs(sum(Q) * 1.0E9) << " [nC]" << endl;
        //os << "* Ekin            =   " << setw(12) << setprecision(5) << eKin_m << " [MeV] dEkin= " << dE_m << " [MeV]" << endl;
        //os << "* rmax            = " << setw(12) << setprecision(5) << rmax_m << " [m]" << endl;
        //os << "* rmin            = " << setw(12) << setprecision(5) << rmin_m << " [m]" << endl;
        //os << "* rms beam size   = " << setw(12) << setprecision(5) << rrms_m << " [m]" << endl;
        //os << "* rms momenta     = " << setw(12) << setprecision(5) << prms_m << " [beta gamma]" << endl;
        //os << "* mean position   = " << setw(12) << setprecision(5) << rmean_m << " [m]" << endl;
        //os << "* mean momenta    = " << setw(12) << setprecision(5) << pmean_m << " [beta gamma]" << endl;
        //os << "* rms emittance   = " << setw(12) << setprecision(5) << eps_m << " (not normalized)" << endl;
        //os << "* rms correlation = " << setw(12) << setprecision(5) << rprms_m << endl;
        //os << "* hr              = " << setw(12) << setprecision(5) << hr_m << " [m]" << endl;
        //os << "* dh              =   " << setw(12) << setprecision(5) << dh_m << " [m]" << endl;
        //os << "* t               =   " << setw(12) << setprecision(5) << getT() << " [s] dT= " << getdT() << " [s]" << endl;
        //os << "* spos            =   " << setw(12) << setprecision(5) << get_sPos() << " [m]" << endl;
        //os << "* ********************************************************************************** " << endl;
        //os.flags(ff);
    }
    return os;
}