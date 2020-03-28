//
// MapAnalyser class, analyses linear maps from OPAL-map.
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.

#ifndef MAP_ANALYSER_H
#define MAP_ANALYSER_H

#include "Classic/FixedAlgebra/FMatrix.h"

#include "Utility/IpplTimings.h"

#include <array>
#include <complex>

/*!
// ------------------------------------------------------------------------
//
// Class: MapAnalyser
//   Organizes the function for a linear map analysis from
//   ThickTracker.
//   Transfer map -> tunes, symplecticity and stability
//   Sigma Matrix -> (not projected) beam emittance
//
// ------------------------------------------------------------------------
//
// $Author: ganz $
//
// ------------------------------------------------------------------------
*/
class MapAnalyser
{
public:
    typedef FMatrix<double, 6, 6> fMatrix_t;
    typedef FMatrix<std::complex<double>, 6, 6> cfMatrix_t;

    MapAnalyser();

    /*!
     * Analyzes the transfer matrix for tunes, symplecticity and stability
     *
     * @param tMatrix Transfer matrix
     */
    void linTAnalyze(const fMatrix_t& /*tMatrix*/) {};


    /*!
     * Analyzes the transfer matrix
     *
     * @ param sigMatrix Sigma Matrix
     */
    void linSigAnalyze(fMatrix_t& /*sigMatrix*/) {};
private:

    /*!
     * Eigen-decomposition of M
     *
     * \f[
     *  \mathbf{M} = \mathbf{E} \mathbf{\Lambda} \mathbf{E}^{-1}
     * \f]
     *
     * @param M Matrix to analyze
     * @param eigenVal EigenValues of \f$\mathbf{M}\f$ \f$\mathbf{\Lambda}\f$
     * @param eigenVec EigenVectors of \f$\mathbf{M}\f$ \f$\mathbf{E}\f$
     * @param invEigenVec inverted EigenVectors of \f$\mathbf{M}\f$ \f$\mathbf{E}^{-1}\f$
     */
    void eigenDecomp_m(const fMatrix_t& M, cfMatrix_t& eigenVal, cfMatrix_t& eigenVec, cfMatrix_t& invEigenVec);

    /*!
     * Returns the block diagonal form rotation matrix \f$ \mathbf{R}\f$. \n
     * This Matrix gets created by applying a symplectic and its inverse Matrix \f$U\f$ on M \f$ \mathbf{R} = \mathbf{UMU}^{-1}\f$
     * @cite Wolski_2005
     *
     * \f[
     *  \mathbf{R}_{\left( \mu_x,\mu_y, \mu_z  \right)}=
     *  \begin{pmatrix}
     *      \mathbf{R}_{2\left( \mu_x\right)} & 0 &0 \\
     *      0 & \mathbf{R}_{2\left( \mu_y\right)} &0 \\
     *      0 & 0 &\mathbf{R}_{2\left( \mu_z\right)} \end{pmatrix}
     *  , \qquad
     *  \mathbf{R_2}_{\left( \alpha \right)}=
     *  \begin{pmatrix}
     *  \cos\left( \alpha \right) & \sin\left( \alpha \right)\\
     *   -\sin\left( \alpha \right) & \cos\left( \alpha \right)\\
     *   \end{pmatrix}
     * \f]
     *
     * @param M Matrix to analyze
     * @param eigenVec EigenVectors of \f$\mathbf{M}\f$
     * @param invEigenVec inverted EigenVectors of \f$\mathbf{M}\f$
     *
     */
    cfMatrix_t getBlockDiagonal_m(const fMatrix_t& M, cfMatrix_t& eigenVecM, cfMatrix_t& invEigenVecM);

    /*!
     * :TODO: WIP Prints phase shift
     *
     * @param Sigma Sigma Matirx
     * @param tM Transfer Matirx
     * @param oldN \f$\mathbf{N}\f$ matrix
     */
    void printPhaseShift_m(fMatrix_t& Sigma, fMatrix_t tM, cfMatrix_t& oldN);

    /*!
    * sets a symplectic \f$ \mathbf{N} \f$ matrix. \n
    * This function is a subfunction of ThickTracker::getBlockDiagonal. \cite Wolski_2005
    *
    * @param M Matrix to analyze
    * @param N symplectic \f$\mathbf{N}\f$ matrix
    * @param invN inverted symplectic \f$\mathbf{N}^{-1}\f$ matrix
    */
    void setNMatrix_m(fMatrix_t& M, cfMatrix_t& N, cfMatrix_t& invN);


    fMatrix_t createRotMatrix_m(std::array<double, 3> phi);

    fMatrix_t createSkewMatrix_m();

    fMatrix_t realPartOfMatrix_m(cfMatrix_t cM);

    fMatrix_t imagPartOfMatrix_m(cfMatrix_t cM);

    cfMatrix_t complexTypeMatrix_m(fMatrix_t M);

    cfMatrix_t invertMatrix_m(const cfMatrix_t &M);

    void rearrangeEigen_m(cfMatrix_t& /*EigenVal*/, cfMatrix_t& /*EigenVec*/) {};

    void normalizeEigen_m(cfMatrix_t& eigenVec, cfMatrix_t& invEigenVec);

private:
    IpplTimings::TimerRef mapAnalysis_m;
    IpplTimings::TimerRef bunchAnalysis_m;
};

#endif
