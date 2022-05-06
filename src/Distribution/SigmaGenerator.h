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

#ifndef SIGMAGENERATOR_H
#define SIGMAGENERATOR_H

#include <array>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "FixedAlgebra/FTps.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Distribution/RealDiracMatrix.h"

class Cyclotron;

class SigmaGenerator
{
public:
    /// Type for storing maps
    typedef RealDiracMatrix::matrix_t matrix_t;
    /// Type for storing the sparse maps
    typedef RealDiracMatrix::sparse_matrix_t sparse_matrix_t;
    /// Type for storing vectors
    typedef RealDiracMatrix::vector_t vector_t;
    /// Container for storing the properties for each angle
    typedef std::vector<double> container_t;
    /// Type of the truncated power series
    typedef FTps<double,2*3> Series;
    /// Type of a map
    typedef FVps<double,2*3> Map;
    /// Type of the Hamiltonian for the cyclotron
    typedef std::function<Series(double,double,double)> Hamiltonian;
    /// Type of the Hamiltonian for the space charge
    typedef std::function<Series(double,double,double)> SpaceCharge;

    /// Constructs an object of type SigmaGenerator
    /*!
     * @param I specifies the current for which a matched distribution should be found, \f$ [I] = A \f$
     * @param ex is the emittance in x-direction (horizontal), \f$ \left[\varepsilon_{x}\right] = \pi\ mm\ mrad  \f$
     * @param ey is the emittance in y-direction (longitudinal), \f$ \left[\varepsilon_{y}\right] = \pi\ mm\ mrad  \f$
     * @param ez is the emittance in z-direction (vertical), \f$ \left[\varepsilon_{z}\right] = \pi\ mm\ mrad  \f$
     * @param E is the energy, \f$ \left[E\right] = MeV \f$
     * @param m is the mass of the particles \f$ \left[m\right] = \frac{MeV}{c^{2}} \f$
     * @param q is the particle charge [e]
     * @param cycl is the cyclotron element
     * @param N is the number of integration steps (closed orbit computation). That's why its also the number
     *    of maps (for each integration step a map)
     * @param Nsectors is the number of sectors that the field map is averaged over (closed orbit computation)
     * @param truncOrder is the truncation order for power series of the Hamiltonian
     * @param write is a boolean (default: true). If true all maps of all iterations are stored, otherwise not.
     */
    SigmaGenerator(double I, double ex, double ey, double ez,
                   double E, double m, double q, const Cyclotron* cycl,
                   unsigned int N, unsigned int Nsectors, unsigned int truncOrder, bool write = true);

    /// Searches for a matched distribution.
    /*!
     * Returns "true" if a matched distribution within the accuracy could be found, returns "false" otherwise.
     * @param accuracy is used for computing the equilibrium orbit (to a certain accuracy)
     * @param maxit is the maximum number of iterations (the program stops if within this value no stationary
     *    distribution was found)
     * @param maxitOrbit is the maximum number of iterations for finding closed orbit
     * @param angle defines the start of the sector (one can choose any angle between 0° and 360°)
     * @param denergy energy step size for closed orbit finder [MeV]
     * @param rguess value of radius for closed orbit finder
     * @param type specifies the magnetic field format (e.g. CARBONCYCL)
     * @param full match over full turn not just single sector
     */
    bool match(double accuracy, unsigned int maxit, unsigned int maxitOrbit,
               Cyclotron* cycl, double denergy, double rguess, bool full);

    /// Block diagonalizes the symplex part of the one turn transfer matrix
    /*! It computes the transformation matrix <b>R</b> and its inverse <b>invR</b>.
     *
     * @param Mturn is a 6x6 dimensional one turn transfer matrix
     * @param R is the 6x6 dimensional transformation matrix (gets computed)
     * @param invR is the 6x6 dimensional inverse transformation (gets computed)
     */
    vector_t decouple(const matrix_t& Mturn, sparse_matrix_t& R, sparse_matrix_t& invR);

    /// Checks if the sigma-matrix is an eigenellipse and returns the L2 error.
    /*!
     * The idea of this function is taken from Dr. Christian Baumgarten's program.
     * @param Mturn is the one turn transfer matrix
     * @param sigma is the sigma matrix to be tested
     */
    double isEigenEllipse(const matrix_t& Mturn, const matrix_t& sigma);

    /// Returns the converged stationary distribution
    matrix_t& getSigma();

    /// Returns the number of iterations needed for convergence (if not converged, it returns zero)
    unsigned int getIterations() const;

    /// Returns the error (if the program didn't converged, one can call this function to check the error)
    double getError() const;

    /*! Returns the emittances (ex,ey,ez) in \f$ \pi\ mm\ mrad \f$ for which
     * the code converged (since the whole simulation is done on normalized emittances)
     */
    std::array<double,3> getEmittances() const;

    const double& getInjectionRadius() const {
        return rinit_m;
    }

    const double& getInjectionMomentum() const {
        return prinit_m;
    }

private:
    /// Stores the value of the current, \f$ \left[I\right] = A \f$
    double I_m;
    /*! Stores the desired emittances,
     * \f$ \left[\varepsilon_{x}\right] = \left[\varepsilon_{y}\right] = \left[\varepsilon_{z}\right] = m \ rad \f$
     */
    std::array<double,3> emittance_m;
    /// Is the orbital frequency, \f$ \left[\omega_{o}\right] = \frac{1}{s} \f$
    double wo_m;
    /// Stores the user-define energy, \f$ \left[E\right] = MeV \f$
    double E_m;
    /*! Relativistic factor (which can be computed out ot the kinetic
     * energy and rest mass (potential energy)), \f$ \left[\gamma\right] = 1 \f$
     */
    double gamma_m;
    /// Relativistic factor squared
    double gamma2_m;
    /// Harmonic number, \f$ \left[N_{h}\right] = 1 \f$
    double nh_m;
    /// Velocity (c/v), \f$ \left[\beta\right] = 1 \f$
    double beta_m;
    /// Is the mass of the particles, \f$ \left[m\right] = \frac{MeV}{c^{2}} \f$
    double mass_m;
    /// Is the particle charge [e]
    double q_m;
    /// Is the number of iterations needed for convergence
    unsigned int niterations_m;
    /// Is true if converged, false otherwise
    bool converged_m;
    /// Minimum energy needed in cyclotron, \f$ \left[E_{min}\right] = MeV \f$
    double Emin_m;
    /// Maximum energy reached in cyclotron, \f$ \left[E_{max}\right] = MeV \f$
    double Emax_m;
    /// Number of integration steps for closed orbit computation
    unsigned int N_m;
    /// Number of (symmetric) sectors the field is averaged over
    unsigned int nSectors_m;
    /// Number of integration steps per sector (--> also: number of maps per sector)
    unsigned int nStepsPerSector_m;

    /// Number of integration steps in total
    unsigned int nSteps_m;

    /// Error of computation
    double error_m;

    /// Truncation order of the power series for the Hamiltonian (cyclotron and space charge)
    unsigned int truncOrder_m;

    /// Decides for writing output (default: true)
    bool write_m;

    /// Stores the stationary distribution (sigma matrix)
    matrix_t sigma_m;

    /// Vector storing the second moment matrix for each angle
    std::vector<matrix_t> sigmas_m;

    /// Stores the Hamiltonian of the cyclotron
    Hamiltonian H_m;

    /// Stores the Hamiltonian for the space charge
    SpaceCharge Hsc_m;

    /// All variables x, px, z, pz, l, delta
    Series x_m, px_m, z_m, pz_m, l_m, delta_m;

    double rinit_m;
    double prinit_m;

    /*! Initializes a first guess of the sigma matrix with the assumption of
     * a spherical symmetric beam (ex = ey = ez). For each angle split the
     * same initial guess is taken.
     *
     * @param nuz is the vertical tune
     * @param ravg is the average radius of the closed orbit
     */
    void initialize(double, double);

    /// Computes the new initial sigma matrix
    /*!
     * @param M is the 6x6 one turn transfer matrix
     * @param R is the transformation matrix
     * @param invR is the inverse transformation matrix
     */
    matrix_t updateInitialSigma(const matrix_t&,
                                   const vector_t&,
                                   sparse_matrix_t&,
                                   sparse_matrix_t&);

    /// Computes new sigma matrices (one for each angle)
    /*!
     * Mscs is a vector of all space charge maps
     * Mcycs is a vector of all cyclotron maps
     */
    void updateSigma(const std::vector<matrix_t>&,
                     const std::vector<matrix_t>&);

    /// Returns the L2-error norm between the old and new sigma-matrix
    /*!
     * @param oldS is the old sigma matrix
     * @param newS is the new sigma matrix
     */
    double L2ErrorNorm(const matrix_t&, const matrix_t&);

    /// Transforms a floating point value to a string
    /*!
     * @param val is the floating point value which is transformed to a string
     */
    std::string float2string(double val);

    /// Called within SigmaGenerator::match().
    /*!
     * @param r_turn is the radius [m]
     * @param peo is the momentum
     * @param h_turn is the inverse bending radius
     * @param fidx_turn is the field index
     * @param ds_turn is the path length element
     */
    void writeOrbitOutput_m(const container_t& r_turn,
                            const container_t& peo,
                            const container_t& h_turn,
                            const container_t& fidx_turn,
                            const container_t& ds_turn);

    void writeMatrix(std::ofstream&, const matrix_t&);

    /// <b>RealDiracMatrix</b>-class member used for decoupling
    RealDiracMatrix rdm_m;


    template<class matrix>
    void reduce(matrix&);

    template<class matrix>
    void expand(matrix&);
};


inline
typename SigmaGenerator::matrix_t&
SigmaGenerator::getSigma()
{
    return sigma_m;
}


inline
unsigned int SigmaGenerator::getIterations() const
{
    return (converged_m) ? niterations_m : 0;
}


inline
double SigmaGenerator::getError() const
{
    return error_m;
}


inline
std::array<double,3> SigmaGenerator::getEmittances() const
{
    double bgam = gamma_m*beta_m;
    return std::array<double,3>{{
            emittance_m[0] / Physics::pi / bgam * Units::m2mm * Units::rad2mrad,
        emittance_m[1] / Physics::pi / bgam * Units::m2mm * Units::rad2mrad,
        emittance_m[2] / Physics::pi / bgam * Units::m2mm * Units::rad2mrad
    }};
}


template<class matrix>
void SigmaGenerator::reduce(matrix& M) {
    /* The 6x6 matrix gets reduced to a 4x4 matrix in the following way:
     *
     * a11 a12 a13 a14 a15 a16
     * a21 a22 a23 a24 a25 a26          a11 a12 a15 a16
     * a31 a32 a33 a34 a35 a36  -->     a21 a22 a25 a26
     * a41 a42 a43 a44 a45 a46          a51 a52 a55 a56
     * a51 a52 a53 a54 a55 a56          a61 a62 a65 a66
     * a61 a62 a63 a64 a65 a66
     */

    // copy x- and l-direction to a 4x4 matrix_t
    matrix_t M4x4(4,4);
    for (unsigned int i = 0; i < 2; ++i) {
        // upper left 2x2 [a11,a12;a21,a22]
        M4x4(i,0) = M(i,0);
        M4x4(i,1) = M(i,1);
        // lower left 2x2 [a51,a52;a61,a62]
        M4x4(i + 2,0) = M(i + 4,0);
        M4x4(i + 2,1) = M(i + 4,1);
        // upper right 2x2 [a15,a16;a25,a26]
        M4x4(i,2) = M(i,4);
        M4x4(i,3) = M(i,5);
        // lower right 2x2 [a55,a56;a65,a66]
        M4x4(i + 2,2) = M(i + 4,4);
        M4x4(i + 2,3) = M(i + 4,5);
    }

    M.resize(4,4,false);
    M.swap(M4x4);
}

template<class matrix>
void SigmaGenerator::expand(matrix& M) {
    /* The 4x4 matrix gets expanded to a 6x6 matrix in the following way:
     *
     *                          a11 a12 0 0 a13 a14
     * a11 a12 a13 a14          a21 a22 0 0 a23 a24
     * a21 a22 a23 a24  -->     0   0   1 0 0   0
     * a31 a32 a33 a34          0   0   0 1 0   0
     * a41 a42 a43 a44          a31 a32 0 0 a33 a34
     *                          a41 a42 0 0 a43 a44
     */

    matrix M6x6 = boost::numeric::ublas::identity_matrix<double>(6,6);

    for (unsigned int i = 0; i < 2; ++i) {
        // upper left 2x2 [a11,a12;a21,a22]
        M6x6(i,0) = M(i,0);
        M6x6(i,1) = M(i,1);
        // lower left 2x2 [a31,a32;a41,a42]
        M6x6(i + 4,0) = M(i + 2,0);
        M6x6(i + 4,1) = M(i + 2,1);
        // upper right 2x2 [a13,a14;a23,a24]
        M6x6(i,4) = M(i,2);
        M6x6(i,5) = M(i,3);
        // lower right 2x2 [a22,a34;a43,a44]
        M6x6(i + 4,4) = M(i + 2,2);
        M6x6(i + 4,5) = M(i + 2,3);
    }

    // exchange
    M.swap(M6x6);
}

#endif