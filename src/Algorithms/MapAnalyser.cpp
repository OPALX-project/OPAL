//
// Source file of the MapAnalyser class, analyses linear maps from OPAL-map.
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.

#include "MapAnalyser.h"

#include <fstream>

#include "Physics/Physics.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>


MapAnalyser::MapAnalyser()
    : mapAnalysis_m(IpplTimings::getTimer("mapAnalysis"))
    , bunchAnalysis_m(IpplTimings::getTimer("bunchAnalysis"))
{ }

//Returns the eigenDecompositon for the fMatrix M = eigenVec * eigenVal * invEigenVec
void MapAnalyser::eigenDecomp_m(const fMatrix_t& M, cfMatrix_t& eigenVal, cfMatrix_t& eigenVec, cfMatrix_t& invEigenVec){

    double data[36];
    int idx, s;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            idx = i * 6 + j;
            data[idx] = M[i][j];
        }
    }

    gsl_matrix_view m = gsl_matrix_view_array(data, 6, 6);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(6);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(6, 6);
    gsl_matrix_complex *eveci = gsl_matrix_complex_alloc(6, 6);
    gsl_permutation * p = gsl_permutation_alloc(6);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(6);

    //get Eigenvalues and Eigenvectors
    gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free(w);
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    for (int i = 0; i < 6; ++i) {
        eigenVal[i][i] = std::complex<double>(
                                              GSL_REAL(gsl_vector_complex_get(eval, i)),
                                              GSL_IMAG(gsl_vector_complex_get(eval, i)));
        for (int j = 0; j < 6; ++j) {
            eigenVec[i][j] = std::complex<double>(
                                                  GSL_REAL(gsl_matrix_complex_get(evec, i, j)),
                                                  GSL_IMAG(gsl_matrix_complex_get(evec, i, j)));
        }
    }

    //invert Eigenvectormatrix
    gsl_linalg_complex_LU_decomp(evec, p, &s);
    gsl_linalg_complex_LU_invert(evec, p, eveci);

    //Create invEigenVecMatrix
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            invEigenVec[i][j] = std::complex<double>(
                                                     GSL_REAL(gsl_matrix_complex_get(eveci, i, j)),
                                                     GSL_IMAG(gsl_matrix_complex_get(eveci, i, j)));
        }
    }

    //free space
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_matrix_complex_free(eveci);
}


//Transforms the Matirx to a block diagonal rotation Matrix
MapAnalyser::cfMatrix_t MapAnalyser::getBlockDiagonal_m(const fMatrix_t& M,
                                                        cfMatrix_t& eigenVecM, cfMatrix_t& invEigenVecM){

    cfMatrix_t cM, qMatrix, invqMatrix, nMatrix, invnMatrix, rMatrix;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            cM[i][j]=std::complex<double>(M[i][j], 0);
        }
    }

    for (int i = 0; i < 6; i = i+2){
        qMatrix[  i][  i] = std::complex<double>(1., 0);
        qMatrix[  i][1+i] = std::complex<double>(0, 1.);
        qMatrix[1+i][  i] = std::complex<double>(1., 0);
        qMatrix[1+i][1+i] = std::complex<double>(0, -1);

        invqMatrix[  i][  i] = std::complex<double>(1., 0);
        invqMatrix[  i][1+i] = std::complex<double>(1., 0);
        invqMatrix[1+i][  i] = std::complex<double>(0., -1.);
        invqMatrix[1+i][1+i] = std::complex<double>(0, 1.);
    }
    qMatrix /= std::sqrt(2.);
    invqMatrix /= std::sqrt(2);

    nMatrix = eigenVecM*qMatrix;
    invnMatrix = invqMatrix* invEigenVecM;


    rMatrix = invnMatrix * cM * nMatrix;

    return rMatrix;
}


void MapAnalyser::printPhaseShift_m(fMatrix_t& Sigma, fMatrix_t tM, cfMatrix_t& oldN){
    cfMatrix_t N1, cinvN, cR, ctM, N2;
    fMatrix_t R1, S, sigmaS;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            ctM[i][j] = std::complex<double>(tM[i][j], 0);
        }
    }

    S = createSkewMatrix_m();
    sigmaS = Sigma*S;

    setNMatrix_m(sigmaS, N2, cinvN);

    std::array<double, 3> phi;

    for (int i = 0; i < 3; i++){
        phi[i] = std::atan(oldN[2*i+1][i].real()/oldN[2*i][2*i].real());
    }

    R1 = createRotMatrix_m(phi);

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            cR[i][j] = std::complex<double>(R1[i][j], 0);
            N1[i][j] = oldN[i][j].real();
        }
    }
}


void MapAnalyser::setNMatrix_m(fMatrix_t& M, cfMatrix_t& N, cfMatrix_t& invN){

    cfMatrix_t eigenValM, eigenVecM, invEigenVecM, eigenVecMT;

    eigenDecomp_m(M, eigenValM, eigenVecM, invEigenVecM);

    cfMatrix_t cM, qMatrix, invqMatrix, nMatrix, invnMatrix, rMatrix;

    //std::ofstream tmap;
    //tmap.open ("TransferMap.txt",std::ios::app);
    //tmap << std::setprecision(16);


    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            cM[i][j] = std::complex<double>(M[i][j], 0);
        }
    }

    for (int i = 0; i < 6; i = i+2){
        qMatrix[  i][  i] = std::complex<double>(1., 0);
        qMatrix[  i][1+i] = std::complex<double>(0, 1.);
        qMatrix[1+i][  i] = std::complex<double>(1., 0);
        qMatrix[1+i][1+i] = std::complex<double>(0, -1);

        invqMatrix[  i][  i] = std::complex<double>(1., 0);
        invqMatrix[  i][1+i] = std::complex<double>(1., 0);
        invqMatrix[1+i][  i] = std::complex<double>(0., -1.);
        invqMatrix[1+i][1+i] = std::complex<double>(0, 1.);
    }
    qMatrix    /= std::sqrt(2);
    invqMatrix /= std::sqrt(2);


    N = eigenVecM*qMatrix;
    invN =  invqMatrix* invEigenVecM;
}


MapAnalyser::fMatrix_t MapAnalyser::createRotMatrix_m(std::array<double, 3> phi){
    fMatrix_t R;

    for (int i = 0; i < 3; i++){
        R[2*i][2*i] = std::cos(phi[1]);
        R[2*i+1][2*i+1] = R[2*i][2*i];
        R[2*i][2*i+1] = std::sin(phi[1]);
        R[2*i+1][2*i] = -R[2*i][2*i+1];
    }
    return R;
}


MapAnalyser::fMatrix_t MapAnalyser::createSkewMatrix_m(){
    fMatrix_t S;

    for (int i = 0; i < 3; i++){
        S[2*i][2*i+1] = 1;
        S[2*i+1][2*i] = -1;
    }
    return S;
}


MapAnalyser::fMatrix_t MapAnalyser::realPartOfMatrix_m(cfMatrix_t cM){
    fMatrix_t M;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            M[i][j] = cM[i][j].real();
        }
    }
    return M;
}

MapAnalyser::fMatrix_t MapAnalyser::imagPartOfMatrix_m(cfMatrix_t cM){
    fMatrix_t M;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            M[i][j] = cM[i][j].imag();
        }
    }
    return M;
}

MapAnalyser::cfMatrix_t MapAnalyser::complexTypeMatrix_m(fMatrix_t M){
    cfMatrix_t cM;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            cM[i][j] = std::complex<double>(M[i][j], 0);
        }
    }
    return cM;
}

MapAnalyser::cfMatrix_t MapAnalyser::invertMatrix_m(const cfMatrix_t& M){

    gsl_set_error_handler_off();
    //gsl_vector_complex *m = gsl_vector_complex_alloc(6);
    gsl_matrix_complex *m = gsl_matrix_complex_alloc(6, 6);
    gsl_matrix_complex *invm = gsl_matrix_complex_alloc(6, 6);
    gsl_permutation * p = gsl_permutation_alloc(6);
    gsl_complex temp;
    int s;


    //Create invEigenVecMatrix
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            GSL_SET_COMPLEX(&temp, std::real(M[i][j]), std::imag(M[i][j]));
            gsl_matrix_complex_set( m, i, j, temp);
        }
    }

    //invert Eigenvectormatrix
    int eigenDecompStatus = gsl_linalg_complex_LU_decomp(m, p, &s);
    if (eigenDecompStatus != 0){
        std::cout<< "This actually works" << std::endl;
        //gsl_set_error_handler (NULL);

    }

    int invertStatus = gsl_linalg_complex_LU_invert(m, p, invm);

    if ( invertStatus ) {
        std::cout << "Error" << std::endl;
        std::exit(1);
    }


    if (invertStatus != 0){
        std::cout<< "This actually works" << std::endl;
        //gsl_set_error_handler (NULL);

    }

    cfMatrix_t invM;
    //Create invEigenVecMatrix
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            invM[i][j] = std::complex<double>(
                                              GSL_REAL(gsl_matrix_complex_get(invm, i, j)),
                                              GSL_IMAG(gsl_matrix_complex_get(invm, i, j)));
        }
    }

    //free space
    gsl_matrix_complex_free(m);
    gsl_matrix_complex_free(invm);
    gsl_permutation_free(p);


    return invM;
}

void MapAnalyser::normalizeEigen_m(cfMatrix_t& eigenVecM, cfMatrix_t& invEigenVecM) {
    //normalize eigen Vectors
    for (int i = 0; i < 6; i++){
        double temp = 0;

        for (int j = 0; j < 6; j += 2){
            temp += 2 * (eigenVecM[j][i] * std::conj(eigenVecM[j+1][i])).imag();
        }
        temp = std::fabs(temp);

        if (temp > 1e-10){
            for (int j = 0; j < 6; j++){
                eigenVecM[j][i] /= std::sqrt(temp);
                invEigenVecM[j][i] /= std::sqrt(temp);
            }
        }
    }
}