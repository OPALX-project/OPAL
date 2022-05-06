//
// Class: SigmaGenerator
// The SigmaGenerator class uses the class <b>ClosedOrbitFinder</b> to get the parameters(inverse bending radius, path length
// field index and tunes) to initialize the sigma matrix.
// The main function of this class is <b>match(double, unsigned int)</b>, where it iteratively tries to find a matched
// distribution for given
// emittances, energy and current. The computation stops when the L2-norm is smaller than a user-defined tolerance. \n
// In default mode it prints all space charge maps, cyclotron maps and second moment matrices. The orbit properties, i.e.
// tunes, average radius, orbit radius, inverse bending radius, path length, field index and frequency error, are printed
// as well.
//
// Copyright (c) 2014, 2018, Matthias Frey, Cristopher Cortes, ETH Zürich
// All rights reserved
//
// Implemented as part of the Semester thesis by Matthias Frey
// "Matched Distributions in Cyclotrons"
//
// Some adaptations done as part of the Bachelor thesis by Cristopher Cortes
// "Limitations of a linear transfer map method for finding matched distributions in high intensity cyclotrons"
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
#include "SigmaGenerator.h"

#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Cyclotron.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include "matrix_vector_operation.h"
#include "ClosedOrbitFinder.h"
#include "MapGenerator.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <utility>

#include <boost/filesystem.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

extern Inform *gmsg;

SigmaGenerator::SigmaGenerator(double I,
                               double ex,
                               double ey,
                               double ez,
                               double E,
                               double m,
                               double q,
                               const Cyclotron* cycl,
                               unsigned int N,
                               unsigned int Nsectors,
                               unsigned int truncOrder,
                               bool write)
    : I_m(I)
    , wo_m(cycl->getRfFrequ()[0]*Units::MHz2Hz/cycl->getCyclHarm()*2.0*Physics::pi)
    , E_m(E)
    , gamma_m(E/m+1.0)
    , gamma2_m(gamma_m*gamma_m)
    , nh_m(cycl->getCyclHarm())
    , beta_m(std::sqrt(1.0-1.0/gamma2_m))
    , mass_m(m)
    , q_m(q)
    , niterations_m(0)
    , converged_m(false)
    , Emin_m(cycl->getFMLowE())
    , Emax_m(cycl->getFMHighE())
    , N_m(N)
    , nSectors_m(Nsectors)
    , nStepsPerSector_m(N/cycl->getSymmetry())
    , nSteps_m(N)
    , error_m(std::numeric_limits<double>::max())
    , truncOrder_m(truncOrder)
    , write_m(write)
    , sigmas_m(nStepsPerSector_m)
    , rinit_m(0.0)
    , prinit_m(0.0)
{
    // minimum beta*gamma
    double bgam = Util::getBetaGamma(Emin_m, mass_m);

    // set emittances (initialization like that due to old compiler version)
    // [ex] = [ey] = [ez] = pi*mm*mrad --> [emittance] = m rad
    // normalized emittance (--> multiply with beta*gamma)
    emittance_m[0] = ex * Physics::pi * bgam * Units::mm2m * Units::mrad2rad;
    emittance_m[1] = ey * Physics::pi * bgam * Units::mm2m * Units::mrad2rad;
    emittance_m[2] = ez * Physics::pi * bgam * Units::mm2m * Units::mrad2rad;

    // Define the Hamiltonian
    Series::setGlobalTruncOrder(truncOrder_m);

    // infinitesimal elements
    x_m = Series::makeVariable(0);
    px_m = Series::makeVariable(1);
    z_m = Series::makeVariable(2);
    pz_m = Series::makeVariable(3);
    l_m = Series::makeVariable(4);
    delta_m = Series::makeVariable(5);

    H_m = [&](double h, double kx, double ky) {
        return 0.5*px_m*px_m + 0.5*kx*x_m*x_m - h*x_m*delta_m +
               0.5*pz_m*pz_m + 0.5*ky*z_m*z_m +
               0.5*delta_m*delta_m/gamma2_m;
    };

    Hsc_m = [&](double sigx, double sigy, double sigz) {
        double m = mass_m * Units::MeV2eV;

        // formula (57)
        double lam = Physics::two_pi*Physics::c / (wo_m * nh_m); // wavelength, [lam] = m
        // [K3] = m
        double K3 = 3.0 * std::abs(q_m) * I_m * lam
                  / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * m
                          * Physics::c * beta_m * beta_m * gamma_m * gamma2_m);

        // formula (30), (31)
        // [sigma(0,0)] = m^{2} --> [sx] = [sy] = [sz] = m
        // In the cyclotron community z is the vertical and y the longitudinal
        // direction.
        // x: horizontal
        // y/l: longitudinal
        // z: vertical
        double sx = std::sqrt(std::abs(sigx));
        double sy = std::sqrt(std::abs(sigy));
        double sz = std::sqrt(std::abs(sigz));

        double tmp = sx * sy;                                           // [tmp] = m^{2}

        double f = std::sqrt(tmp) / (3.0 * gamma_m * sz);               // [f] = 1
        double kxy = K3 * std::abs(1.0 - f) / ((sx + sy) * sz); // [kxy] = 1/m

        double Kx = kxy / sx;
        double Ky = kxy / sy;
        double Kz = K3 * f / (tmp * sz);

        return -0.5 * Kx * x_m * x_m
               -0.5 * Kz * z_m * z_m
               -0.5 * Ky * l_m * l_m * gamma2_m;
     };
}


bool SigmaGenerator::match(double accuracy,
                           unsigned int maxit,
                           unsigned int maxitOrbit,
                           Cyclotron* cycl,
                           double denergy,
                           double rguess,
                           bool full)
{
    /* compute the equilibrium orbit for energy E_
     * and get the the following properties:
     * - inverse bending radius h
     * - step sizes of path ds
     * - tune nuz
     */

    try {
        if ( !full )
            nSteps_m = nStepsPerSector_m;

        // object for space charge map and cyclotron map
        MapGenerator<double,
                     unsigned int,
                     Series,
                     Map,
                     Hamiltonian,
                     SpaceCharge> mapgen(nSteps_m);

        // compute cyclotron map and space charge map for each angle and store them into a vector
        std::vector<matrix_t> Mcycs(nSteps_m), Mscs(nSteps_m);

        ClosedOrbitFinder<double, unsigned int,
            boost::numeric::odeint::runge_kutta4<container_t> > cof(mass_m, q_m, N_m, cycl, false, nSectors_m);

        if ( !cof.findOrbit(accuracy, maxitOrbit, E_m, denergy, rguess) ) {
            throw OpalException("SigmaGenerator::match()",
                                "Closed orbit finder didn't converge.");
        }

        cof.computeOrbitProperties(E_m);

        double angle = cycl->getPHIinit();
        container_t h    = cof.getInverseBendingRadius(angle);
        container_t r    = cof.getOrbit(angle);
        container_t peo  = cof.getMomentum(angle);
        container_t ds   = cof.getPathLength(angle);
        container_t fidx = cof.getFieldIndex(angle);

        // write properties to file
        writeOrbitOutput_m(r, peo, h, fidx, ds);

        rinit_m = r[0];
        prinit_m = peo[0];

        std::string fpath = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps"
        });
        if (!boost::filesystem::exists(fpath)) {
            boost::filesystem::create_directory(fpath);
        }

        std::pair<double,double> tunes = cof.getTunes();
        double ravg = cof.getAverageRadius();

        // write to terminal
        *gmsg << "* ----------------------------" << endl
                << "* Closed orbit info:" << endl
                << "*" << endl
                << "* average radius: " << ravg << " [m]" << endl
                << "* initial radius: " << r[0] << " [m]" << endl
                << "* initial momentum: " << peo[0] << " [Beta Gamma]" << endl
                << "* frequency error: " << cof.getFrequencyError()*100 <<" [ % ] "<< endl
                << "* horizontal tune: " << tunes.first << endl
                << "* vertical tune: " << tunes.second << endl
                << "* ----------------------------" << endl << endl;

        // initialize sigma matrices (for each angle one) (first guess)
        initialize(tunes.second,ravg);

        // for writing
        std::ofstream writeMturn, writeMcyc, writeMsc;

        if (write_m) {

            std::string energy = float2string(E_m);

            std::string fname = Util::combineFilePath({
                OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                "maps",
                "OneTurnMapsForEnergy" + energy + "MeV.dat"
            });

            writeMturn.open(fname, std::ios::out);

            fname = Util::combineFilePath({
                OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                "maps",
                "SpaceChargeMapPerAngleForEnergy" + energy + "MeV_iteration_0.dat"
            });

            writeMsc.open(fname, std::ios::out);

            fname = Util::combineFilePath({
                OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                "maps",
                "CyclotronMapPerAngleForEnergy" + energy + "MeV.dat"
            });

            writeMcyc.open(fname, std::ios::out);
        }

        // calculate only for a single sector (a nSector_-th) of the whole cyclotron
        for (unsigned int i = 0; i < nSteps_m; ++i) {
            Mcycs[i] = mapgen.generateMap(H_m(h[i],
                                              h[i]*h[i]+fidx[i],
                                              -fidx[i]),
                                          ds[i],truncOrder_m);


            Mscs[i]  = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),
                                                sigmas_m[i](2,2),
                                                sigmas_m[i](4,4)),
                                          ds[i],truncOrder_m);

            writeMatrix(writeMcyc, Mcycs[i]);
            writeMatrix(writeMsc, Mscs[i]);
        }

        if (write_m)
            writeMsc.close();

        // one turn matrix
        mapgen.combine(Mscs,Mcycs);
        matrix_t Mturn = mapgen.getMap();

        writeMatrix(writeMturn, Mturn);

        // (inverse) transformation matrix
        sparse_matrix_t R(6, 6), invR(6, 6);

        // new initial sigma matrix
        matrix_t newSigma(6,6);

        // for exiting loop
        bool stop = false;

        double weight = 0.5;

        while (error_m > accuracy && !stop) {
            // decouple transfer matrix and compute (inverse) tranformation matrix
           vector_t eigen = decouple(Mturn, R,invR);

            // construct new initial sigma-matrix
            newSigma = updateInitialSigma(Mturn, eigen, R, invR);

            // compute new sigma matrices for all angles (except for initial sigma)
            updateSigma(Mscs,Mcycs);

            // compute error with mm^2 and (mrad)^2
            error_m = L2ErrorNorm(sigmas_m[0] * 1e6, newSigma * 1e6);

            // write new initial sigma-matrix into vector
            sigmas_m[0] = weight*newSigma + (1.0-weight)*sigmas_m[0];

            if (write_m) {

                std::string energy = float2string(E_m);

                std::string fname = Util::combineFilePath({
                        OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                        "maps",
                        "SpaceChargeMapPerAngleForEnergy" + energy + "MeV_iteration_"
                        + std::to_string(niterations_m + 1)
                        + ".dat"
                });

                writeMsc.open(fname, std::ios::out);
            }

            // compute new space charge maps
            for (unsigned int i = 0; i < nSteps_m; ++i) {
                Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),
                                                    sigmas_m[i](2,2),
                                                    sigmas_m[i](4,4)),
                                                ds[i],truncOrder_m);

                writeMatrix(writeMsc, Mscs[i]);
            }

            if (write_m) {
                writeMsc.close();
            }

            // construct new one turn transfer matrix M
            mapgen.combine(Mscs,Mcycs);
            Mturn = mapgen.getMap();

            writeMatrix(writeMturn, Mturn);

            // check if number of iterations has maxit exceeded.
            stop = (niterations_m++ > maxit);
        }

        // store converged sigma-matrix
        sigma_m.resize(6,6,false);
        sigma_m.swap(newSigma);

        // returns if the sigma matrix has converged
        converged_m = error_m < accuracy;

        // Close files
        if (write_m) {
            writeMturn.close();
            writeMcyc.close();
        }

    } catch(const std::exception& e) {
        throw OpalException("SigmaGenerator::match()", e.what());
    }

    if ( converged_m && write_m ) {
        // write tunes
        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "MatchedDistributions.dat"
        });

        std::ofstream writeSigmaMatched(fname, std::ios::out);

        std::array<double,3> emit = this->getEmittances();

        writeSigmaMatched << "* Converged (Ex, Ey, Ez) = (" << emit[0]
                          << ", " << emit[1] << ", " << emit[2]
                          << ") pi mm mrad for E= " << E_m << " (MeV)"
                          << std::endl << "* Sigma-Matrix " << std::endl;

        for(unsigned int i = 0; i < sigma_m.size1(); ++ i) {
            writeSigmaMatched << std::setprecision(4)  << std::setw(11)
                              << sigma_m(i,0);
            for(unsigned int j = 1; j < sigma_m.size2(); ++ j) {
                writeSigmaMatched << " & " <<  std::setprecision(4)
                                  << std::setw(11) << sigma_m(i,j);
            }
            writeSigmaMatched << " \\\\" << std::endl;
        }

        writeSigmaMatched << std::endl;
        writeSigmaMatched.close();
    }

    return converged_m;
}


typename SigmaGenerator::vector_t
SigmaGenerator::decouple(const matrix_t& Mturn,
                         sparse_matrix_t& R,
                         sparse_matrix_t& invR)
{
    // copy one turn matrix
    matrix_t M(Mturn);

    // reduce 6x6 matrix to 4x4 matrix
    reduce<matrix_t>(M);

    // compute symplex part
    matrix_t Ms = rdm_m.symplex(M);

    // diagonalize and compute transformation matrices
    rdm_m.diagonalize(Ms,R,invR);

    /*
     * formula (38) in paper of Dr. Christian Baumgarten:
     * Geometrical method of decoupling
     *
     *          [0,     alpha,  0,      0;
     * F_{d} =  -beta,  0,      0,      0;
     *          0,      0,      0,      gamma;
     *          0,      0,      -delta, 0]
     *
     *
     */
    vector_t eigen(4);
    eigen(0) =   Ms(0,1);       // alpha
    eigen(1) = - Ms(1,0);       // beta
    eigen(2) =   Ms(2,3);       // gamma
    eigen(3) = - Ms(3,2);       // delta
    return eigen;
}


double SigmaGenerator::isEigenEllipse(const matrix_t& Mturn,
                                      const matrix_t& sigma)
{
    // compute sigma matrix after one turn
    matrix_t newSigma = matt_boost::gemmm<matrix_t>(Mturn,
                                                    sigma,
                                                    boost::numeric::ublas::trans(Mturn));

    // return L2 error
    return L2ErrorNorm(sigma,newSigma);
}


void SigmaGenerator::initialize(double nuz, double ravg)
{
    /*
     * The initialization is based on the analytical solution of
     * using a spherical symmetric beam in the paper:
     * Transverse-longitudinal coupling by space charge in cyclotrons
     * by Dr. Christian Baumgarten
     * (formulas: (46), (56), (57))
     */


    /* Units:
     * ----------------------------------------------
     * [wo] = 1/s
     * [nh] = 1
     * [q0] = 1 e
     * [I] = A
     * [eps0] = (A*s)^{2}/(N*m^{2})
     * [E0] = MeV/(c^{2}) (with speed of light c)
     * [beta] = 1
     * [gamma] = 1
     * [m] = eV/c^2
     *
     * [lam] = m
     * [K3] = m
     * [alpha] = 1
     * ----------------------------------------------
     */

    // helper constants
    double invbg = 1.0 / (beta_m * gamma_m);

    double mass = mass_m * Units::MeV2eV;

    // emittance [ex] = [ey] = [ez] = m rad
    double ex = emittance_m[0] * invbg;                        // [ex] = m rad
    double ey = emittance_m[1] * invbg;                        // [ey] = m rad
    double ez = emittance_m[2] * invbg;                        // [ez] = m rad

    // initial guess of emittance, [e] = m rad
    double guessedEmittance = std::cbrt(ex * ey * ez);             // cbrt computes cubic root (C++11) <cmath>

    // cyclotron radius [rcyc] = m
    double rcyc = ravg / beta_m;

    // "average" inverse bending radius
    double avgInverseBendingRadius = 1.0 / ravg;

    // formula (57)
    double lam = Physics::two_pi * Physics::c / (wo_m * nh_m); // wavelength, [lam] = m

    // m * c^3 --> c^2 in [m] = eV / c^2 cancel --> m * c in denominator
    double K3 = 3.0 *  std::abs(q_m) * I_m * lam
              / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * mass
                      * Physics::c * beta_m * beta_m * gamma2_m * gamma_m);    // [K3] = m

    // c in denominator cancels with unit of [m] = eV / c^2 --> we need to multiply
    // with c in order to get dimensionless quantity
    double alpha =  (std::abs(q_m) * Physics::mu_0 * I_m * Physics::c
                     / (5.0 * std::sqrt(10.0) * mass * gamma_m * nh_m)
                     * std::sqrt(rcyc * rcyc * rcyc / (std::pow(guessedEmittance, 3))));   // [alpha] = 1/rad --> [alpha] = 1

    double sig0 = std::sqrt(2.0 * rcyc * guessedEmittance) / gamma_m;                     // [sig0] = m*sqrt(rad) --> [sig0] = m

    // formula (56)
    double sig;                                     // [sig] = m
    if (alpha >= 2.5) {
        sig = sig0 * std::cbrt(1.0 + alpha);            // cbrt computes cubic root (C++11) <cmath>
    } else if (alpha >= 0) {
        sig = sig0 * (1 + alpha * (0.25 - 0.03125 * alpha));
    } else {
        throw OpalException("SigmaGenerator::initialize()",
                            "Negative alpha value: " + std::to_string(alpha) + " < 0");
    }

    // K = Kx = Ky = Kz
    double K = K3 * gamma_m / (3.0 * sig * sig * sig);   // formula (46), [K] = 1/m^{2}
    double kx = std::pow(avgInverseBendingRadius, 2) * gamma2_m;// formula (46) (assumption of an isochronous cyclotron), [kx] = 1/m^{2}

    double a = 0.5 * kx - K;    // formula (22) (with K = Kx = Kz), [a] = 1/m^{2}
    double b = K * K;           // formula (22) (with K = Kx = Kz and kx = h^2*gamma^2), [b] = 1/m^{4}


    // b must be positive, otherwise no real-valued frequency
    if (b < 0)
        throw OpalException("SigmaGenerator::initialize()",
                            "Negative value --> No real-valued frequency.");

    double tmp = a * a - b;           // [tmp] = 1/m^{4}
    if (tmp < 0)
        throw OpalException("SigmaGenerator::initialize()",
                            "Square root of negative number.");

    tmp = std::sqrt(tmp);               // [tmp] = 1/m^{2}

    if (a < tmp)
        throw OpalException("SigmaGenerator::initialize()",
                            "Square root of negative number.");

    if (std::pow(avgInverseBendingRadius, 2) * nuz * nuz <= K)
        throw OpalException("SigmaGenerator::initialize()",
                            "h^{2} * nu_{z}^{2} <= K (i.e. square root of negative number)");

    double Omega = std::sqrt(a + tmp);                        // formula (22), [Omega] = 1/m
    double omega = std::sqrt(a - tmp);                        // formula (22), [omega] = 1/m

    double A = avgInverseBendingRadius / (Omega * Omega + K); // formula (26), [A] = m
    double B = avgInverseBendingRadius / (omega * omega + K); // formula (26), [B] = m
    double invAB = 1.0 / (B - A);                 // [invAB] = 1/m

    // construct initial sigma-matrix (formula (29, 30, 31)
    matrix_t sigma = boost::numeric::ublas::zero_matrix<double>(6);

    // formula (30), [sigma(0,0)] = m^2 rad
    sigma(0,0) = invAB * (B * ex / Omega + A * ez / omega);

    // [sigma(0,5)] = [sigma(5,0)] = m rad
    sigma(0,5) = sigma(5,0) = invAB * (ex / Omega + ez / omega);

    // [sigma(1,1)] = rad
    sigma(1,1) = invAB * (B * ex * Omega + A * ez * omega);

    // [sigma(1,4)] = [sigma(4,1)] = m rad
    sigma(1,4) = sigma(4,1) = invAB * (ex * Omega+ez * omega) / (K * gamma2_m);

    // formula (31), [sigma(2,2)] = m rad
    sigma(2,2) = ey / (std::sqrt(std::pow(avgInverseBendingRadius,2) * nuz * nuz - K));

    sigma(3,3) = (ey * ey) / sigma(2,2);

    // [sigma(4,4)] = m^{2} rad
    sigma(4,4) = invAB * (A * ex * Omega + B * ez * omega) / (K * gamma2_m);

    // formula (30), [sigma(5,5)] = rad
    sigma(5,5) = invAB * (ex / (B * Omega) + ez / (A * omega));

    // fill in initial guess of the sigma matrix (for each angle the same guess)
    sigmas_m.resize(nSteps_m);
    for (typename std::vector<matrix_t>::iterator it = sigmas_m.begin(); it != sigmas_m.end(); ++it) {
        *it = sigma;
    }

    if (write_m) {
        std::string energy = float2string(E_m);

        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps",
            "InitialSigmaPerAngleForEnergy" + energy + "MeV.dat"
        });

        std::ofstream writeInit(fname, std::ios::out);
        writeMatrix(writeInit, sigma);
        writeInit.close();
    }
}


typename SigmaGenerator::matrix_t
SigmaGenerator::updateInitialSigma(const matrix_t& M,
                                   const vector_t& eigen,
                                   sparse_matrix_t& R,
                                   sparse_matrix_t& invR)
{
    /*
     * This function is based on the paper of Dr. Christian Baumgarten:
     * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons (2012)
     */

    /*
     * Function input:
     * - M: one turn transfer matrix
     * - eigen = {alpha, beta, gamma, delta}
     * - R: transformation matrix (in paper: E)
     * - invR: inverse transformation matrix (in paper: E^{-1}
     */

    // normalize emittances
    double invbg = 1.0 / (beta_m * gamma_m);
    double ex = emittance_m[0] * invbg;
    double ey = emittance_m[1] * invbg;
    double ez = emittance_m[2] * invbg;

    // alpha^2-beta*gamma = 1

    /* 0        eigen(0) 0        0
     * eigen(1) 0        0        0
     * 0        0        0        eigen(2)
     * 0        0        eigen(3) 0
     *
     * M = cos(mux)*[1, 0; 0, 1] + sin(mux)*[alpha, beta; -gamma, -alpha], Book, p. 242
     *
     * -----------------------------------------------------------------------------------
     * X-DIRECTION and L-DIRECTION
     * -----------------------------------------------------------------------------------
     * --> eigen(0) = sin(mux)*betax
     * --> eigen(1) = -gammax*sin(mux)
     *
     * What is sin(mux)?   --> alphax = 0 --> -alphax^2+betax*gammax = betax*gammax = 1
     *
     * Convention: betax > 0
     *
     * 1) betax = 1/gammax
     * 2) eig0 = sin(mux)*betax
     * 3) eig1 = -gammax*sin(mux)
     *
     * eig0 = sin(mux)/gammax
     * eig1 = -gammax*sin(mux) <--> 1/gammax = -sin(mux)/eig1
     *
     * eig0 = -sin(mux)^2/eig1 --> -sin(mux)^2 = eig0*eig1      --> sin(mux) = sqrt(-eig0*eig1)
     *                                                          --> gammax = -eig1/sin(mux)
     *                                                          --> betax = eig0/sin(mux)
     */


    // x-direction
    //double alphax = 0.0;
    double betax  = std::sqrt(std::fabs(eigen(0) / eigen(1)));
    double gammax = 1.0 / betax;

    // l-direction
    //double alphal = 0.0;
    double betal  = std::sqrt(std::fabs(eigen(2) / eigen(3)));
    double gammal = 1.0 / betal;

    /*
     * -----------------------------------------------------------------------------------
     * Z-DIRECTION
     * -----------------------------------------------------------------------------------
     *
     * m22 m23
     * m32 m33
     *
     * m22 = cos(muz) + alpha*sin(muz)
     * m33 = cos(muz) - alpha*sin(muz)
     *
     * --> cos(muz) = 0.5*(m22 + m33)
     *     sin(muz) = sign(m32)*sqrt(1-cos(muz)^2)
     *
     * m22-m33 = 2*alpha*sin(muz) --> alpha = 0.5*(m22-m33)/sin(muz)
     *
     * m23 = betaz*sin(muz)     --> betaz = m23/sin(muz)
     * m32 = -gammaz*sin(muz)   --> gammaz = -m32/sin(muz)
     */

    double cosz = 0.5 * (M(2,2) + M(3,3));

    double sign = (std::signbit(M(2,3))) ? double(-1) : double(1);

    double invsinz = sign / std::sqrt(std::fabs( 1.0 - cosz * cosz));

    double alphaz = 0.5 * (M(2,2) - M(3,3)) * invsinz;
    double betaz  =   M(2,3) * invsinz;
    double gammaz = - M(3,2) * invsinz;

    // Convention beta>0
    if (std::signbit(betaz))    // singbit = true if beta<0, else false
        betaz  *= -1.0;

    // diagonal matrix with eigenvalues
    matrix_t D = boost::numeric::ublas::zero_matrix<double>(6,6);
    // x-direction
    D(0,1) =   betax  * ex;
    D(1,0) = - gammax * ex;
    // z-direction
    D(2,2) =   alphaz * ey;
    D(3,3) = - alphaz * ey;
    D(2,3) =   betaz  * ey;
    D(3,2) = - gammaz * ey;
    // l-direction
    D(4,5) =   betal  * ez;
    D(5,4) = - gammal * ez;

    // expand 4x4 transformation matrices to 6x6
    expand<sparse_matrix_t>(R);
    expand<sparse_matrix_t>(invR);

    // symplectic matrix
    sparse_matrix_t S(6,6,6);
    S(0,1) = S(2,3) = S(4,5) = 1;
    S(1,0) = S(3,2) = S(5,4) = -1;

    matrix_t sigma = matt_boost::gemmm<matrix_t>(-invR,D,R);
    sigma = boost::numeric::ublas::prod(sigma,S);

    if (write_m) {
        std::string energy = float2string(E_m);

        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps",
            "InitialSigmaPerAngleForEnergy" + energy + "MeV.dat"
        });

        std::ofstream writeInit(fname, std::ios::app);
        writeMatrix(writeInit, sigma);
        writeInit.close();
    }

    return sigma;
}


void SigmaGenerator::updateSigma(const std::vector<matrix_t>& Mscs,
                                 const std::vector<matrix_t>& Mcycs)
{
    std::ofstream writeSigma;

    if (write_m) {
        std::string energy = float2string(E_m);

        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps",
            "SigmaPerAngleForEnergy" + energy + "MeV_iteration_"
            + std::to_string(niterations_m + 1)
            + ".dat"
        });

        writeSigma.open(fname,std::ios::out);
    }

    // initial sigma is already computed
    writeMatrix(writeSigma, sigmas_m[0]);

    for (unsigned int i = 1; i < nSteps_m; ++i) {
        // transfer matrix for one angle
        matrix_t M = boost::numeric::ublas::prod(Mscs[i - 1],Mcycs[i - 1]);
        // transfer the matrix sigma
        sigmas_m[i] = matt_boost::gemmm<matrix_t>(M,sigmas_m[i - 1],
                                                     boost::numeric::ublas::trans(M));

        writeMatrix(writeSigma, sigmas_m[i]);
    }

    if (write_m) {
        writeSigma.close();
    }
}


double SigmaGenerator::L2ErrorNorm(const matrix_t& oldS,
                                   const matrix_t& newS)
{
    // compute difference
    matrix_t diff = newS - oldS;

    // sum squared error up and take square root
    return std::sqrt(std::inner_product(diff.data().begin(),
                                        diff.data().end(),
                                        diff.data().begin(),
                                        0.0));
}


std::string SigmaGenerator::float2string(double val) {
    std::ostringstream out;
    out << std::setprecision(6) << val;
    return out.str();
}


void SigmaGenerator::writeOrbitOutput_m(
    const container_t& r,
    const container_t& peo,
    const container_t& h,
    const container_t& fidx,
    const container_t& ds)
{
    if (!write_m)
        return;

    std::string energy = float2string(E_m);
    std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "PropertiesForEnergy" + energy + "MeV.dat"
    });

    std::ofstream writeProperties(fname, std::ios::out);

    writeProperties << std::left
                    << std::setw(25) << "orbit radius"
                    << std::setw(25) << "orbit momentum"
                    << std::setw(25) << "inverse bending radius"
                    << std::setw(25) << "field index"
                    << std::setw(25) << "path length" << std::endl;

    for (unsigned int i = 0; i < r.size(); ++i) {
        writeProperties << std::setprecision(10) << std::left
                        << std::setw(25) << r[i]
                        << std::setw(25) << peo[i]
                        << std::setw(25) << h[i]
                        << std::setw(25) << fidx[i]
                        << std::setw(25) << ds[i] << std::endl;
    }
    writeProperties.close();
}


void SigmaGenerator::writeMatrix(std::ofstream& os, const matrix_t& m) {
    if (!write_m)
        return;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j)
            os << m(i, j) << " ";
    }
    os << std::endl;
}
