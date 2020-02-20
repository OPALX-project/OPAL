// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
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

//Analyzes a TransferMap for the tunes, symplecticity and stability
void MapAnalyser::linTAnalyze(const fMatrix_t& tMatrix){
    std::ofstream tplot;
    tplot.open ("tunePlot.txt",std::ios::app);
    tplot << std::setprecision(16);
    //================================

    IpplTimings::startTimer(mapAnalysis_m);


    fMatrix_t blocktMatrix, scalingMatrix, invScalingMatrix;

    //get the EigenValues/-Vectors
    cfMatrix_t eigenValM, eigenVecM, invEigenVecM;
    eigenDecomp_m(tMatrix, eigenValM, eigenVecM, invEigenVecM);

    //normalize eigen Vectors
    for (int i=0; i <6 ; i++){
        double temp=0;

        for (int j=0; j <6 ; j+=2){
            temp += 2*(eigenVecM[j][i] * std::conj(eigenVecM[j+1][i])).imag();
        }
        temp=std::fabs(temp);

        if (temp > 1e-10){
            //TODO is it neccessary to normalize the eigenValues?
            //eigenValM[i][i] *= std::sqrt(temp);
            for (int j=0; j <6 ; j++){
                eigenVecM[j][i] /= std::sqrt(temp);
                invEigenVecM[j][i] /= std::sqrt(temp);
            }
        }
    }

    //TODO pack transposition in a neat function
    cfMatrix_t eigenVecMT;
    for (int i=0; i <6 ; i++){
        for (int j=0; j <6 ; j++){
            eigenVecMT[i][j]= eigenVecM[j][i];
        }
    }

    double prec =0.01;

    int idx=0;
    std::array<int,6> pairs;

    fMatrix_t skewMatrix = createSkewMatrix_m();
    cfMatrix_t cSkewMatrix = complexTypeMatrix_m(skewMatrix);

    fMatrix_t K = imagPartOfMatrix_m(eigenVecMT *cSkewMatrix* eigenVecM);

    for (int i=0; i<6;i++) pairs[i]=-1;

    //Printer
    for (int i=0; i <6 ; i++){
        for (int j=0; j <6 ; j++){
            if (K[i][j]>1-prec && K[i][j]<1+prec){
                pairs[idx]=i;
                pairs[idx+1]=j;
                idx+=2;

                if (idx==6) break;
            }
        }
    }


    //Fill empty elements in pairs
    std::array<int,6>::iterator pairsIt;

    //Compare eigenvalues
    int negidx=6-1;
    for (int i=0; i<6; i++) {
        pairsIt = std::find (pairs.begin(), pairs.end(), i);
        if (pairsIt != pairs.end()){
            continue;
        }

        for (int j=0; i!=j && j <6 ; j++){
            double diff = std::abs( eigenValM[i][i] - eigenValM[j][j]);
            if (diff< 0.001){
                pairs[negidx]=i;
                pairs[negidx-1]=j;
                negidx-=2;
            }
        }
    }

    //Fill the paris vector
    for (int i=0; i < 6 && idx < 6; i++){
        pairsIt = std::find (pairs.begin(), pairs.end(), i);

        if (pairsIt != pairs.end()){
            continue;
        } else{
            pairs[idx]=i;
            idx++;
        }
    }

    // :FIXME: why commented out? To be removed?
    cfMatrix_t tempM = eigenVecM ;
    //cfMatrix_t tempInvM = eigenVecM ;
    cfMatrix_t tempValM = eigenValM ;


    for (int i=0; i <6 ; i++){
        eigenValM[i][i]=tempValM[pairs[i]][pairs[i]];
        for (int j=0; j <6 ; j++){
            eigenVecM[j][i]=tempM[j][pairs[i]];
            //eigenVecMT[i][j]=eigenVecM[j][i];
            //invEigenVecM[i][j]=tempInvM[pairs[i]][j];
        }
    }

    invEigenVecM= invertMatrix_m(eigenVecM);
    cfMatrix_t cblocktMatrix = getBlockDiagonal_m(tMatrix, eigenVecM, invEigenVecM);

    FVector<std::complex<double>, 3> betaTunes, betaTunes2;
    FVector<double, 3> betaTunes3;

    // :FIXME: why commented out
    //rearrangeEigen(eigenValM, eigenVecM);

    for (int i = 0; i < 3; i++){
        betaTunes[i]=std::log(eigenValM[i*2][i*2])/ (2*Physics::pi * std::complex<double>(0, 1));
        betaTunes[i]= std::asin(cblocktMatrix[i*2][i*2+1])/(std::complex<double>(2*Physics::pi, 0));

        betaTunes2[i]= std::acos(cblocktMatrix[i*2][i*2])/(std::complex<double>(2*Physics::pi, 0));
        betaTunes2[i]= std::acos(cblocktMatrix[i*2][i*2].real())/(std::complex<double>(2*Physics::pi, 0));

        betaTunes3[i]= std::acos(eigenValM[i*2][i*2].real()) / (2*Physics::pi);

        double lenEigenVec = 0;
        for (int j = 0; j < 3; j++){
            lenEigenVec += std::abs(eigenVecM[i][j]);

        }
        lenEigenVec = std::sqrt(lenEigenVec);
        tplot<<"1: "<<betaTunes[i] <<std::endl;
        tplot<<"2: "<<betaTunes2[i]<< std::endl;
        //tplot<<"3: "<<lenEigenVec<< std::endl;
        tplot<<"3: ("<<betaTunes3[i]<< ",0)"<< std::endl;

        IpplTimings::stopTimer(mapAnalysis_m);
        //TODO: do something with the tunes etc
        //:FIXME: are there some plans what to do?
    }

}


//Analyses the Sigma Matrix
void MapAnalyser::linSigAnalyze(fMatrix_t& sigMatrix){
    IpplTimings::startTimer(bunchAnalysis_m);
    fMatrix_t sigmaBlockM;

    fMatrix_t skewMatrix;
    for (int i=0; i <6 ; i=i+2){
        skewMatrix[0+i][1+i]=1.;
        skewMatrix[1+i][0+i]=-1.;
    }

    sigMatrix=sigMatrix*skewMatrix;

    cfMatrix_t eigenValM, eigenVecM, invEigenVecM;

    eigenDecomp_m(sigMatrix,eigenValM, eigenVecM, invEigenVecM);

    cfMatrix_t cblocksigM = getBlockDiagonal_m(sigMatrix, eigenVecM,invEigenVecM);


    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            sigmaBlockM[i][j]=cblocksigM[i][j].real();
        }

    }

    cfMatrix_t teigenValM, teigenVecM, tinvEigenVecM;
    eigenDecomp_m(sigmaBlockM, teigenValM, teigenVecM, tinvEigenVecM);

    for (int i=0; i<6; i++){

    }


    //TODO: Where to go with EigenValues?
    IpplTimings::stopTimer(bunchAnalysis_m);
}


//Returns the eigenDecomositon for the fMatrix M = eigenVec * eigenVal * invEigenVec
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

    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            cM[i][j]=std::complex<double>(M[i][j],0);
        }
    }

    for (int i=0; i <6 ; i=i+2){
        qMatrix[0+i][0+i]=std::complex<double>(1.,0);
        qMatrix[0+i][1+i]=std::complex<double>(0,1.);
        qMatrix[1+i][0+i]=std::complex<double>(1.,0);
        qMatrix[1+i][1+i]=std::complex<double>(0,-1);

        invqMatrix[0+i][0+i]=std::complex<double>(1.,0);
        invqMatrix[0+i][1+i]=std::complex<double>(1.,0);
        invqMatrix[1+i][0+i]=std::complex<double>(0.,-1.);
        invqMatrix[1+i][1+i]=std::complex<double>(0,1.);
    }
    qMatrix/=std::sqrt(2.);
    invqMatrix/=std::sqrt(2);

    nMatrix=eigenVecM*qMatrix;
    invnMatrix= invqMatrix* invEigenVecM;


    rMatrix= invnMatrix * cM * nMatrix;

    return rMatrix;

}


void MapAnalyser::printPhaseShift_m(fMatrix_t& Sigma, fMatrix_t tM, cfMatrix_t& oldN){
    cfMatrix_t N1, cinvN, cR, ctM, N2;
    fMatrix_t R1, S, sigmaS;

    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            ctM[i][j]=std::complex<double>(tM[i][j],0);
        }
    }

    S = createSkewMatrix_m();
    sigmaS = Sigma*S;

    setNMatrix_m(sigmaS, N2, cinvN);

    std::array<double, 3> phi;

    for (int i = 0; i < 3; i++){
        phi[i] = std::atan(oldN[2*i+1][i].real()/oldN[2*i][2*i].real());
    }


    R1=createRotMatrix_m(phi);

    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            cR[i][j]=std::complex<double>(R1[i][j],0);
            N1[i][j]=oldN[i][j].real();
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


    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            cM[i][j]=std::complex<double>(M[i][j],0);
        }
    }

    for (int i=0; i <6 ; i=i+2){
        qMatrix[0+i][0+i]=std::complex<double>(1.,0);
        qMatrix[0+i][1+i]=std::complex<double>(0,1.);
        qMatrix[1+i][0+i]=std::complex<double>(1.,0);
        qMatrix[1+i][1+i]=std::complex<double>(0,-1);

        invqMatrix[0+i][0+i]=std::complex<double>(1.,0);
        invqMatrix[0+i][1+i]=std::complex<double>(1.,0);
        invqMatrix[1+i][0+i]=std::complex<double>(0.,-1.);
        invqMatrix[1+i][1+i]=std::complex<double>(0,1.);
    }
    qMatrix/=std::sqrt(2.);
    invqMatrix/=std::sqrt(2);


    N=eigenVecM*qMatrix;
    invN= invqMatrix* invEigenVecM;
}


MapAnalyser::fMatrix_t MapAnalyser::createRotMatrix_m(std::array<double,3> phi){
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
    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            M[i][j]=cM[i][j].real();
        }
    }
    return M;
}

MapAnalyser::fMatrix_t MapAnalyser::imagPartOfMatrix_m(cfMatrix_t cM){
    fMatrix_t M;
    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            M[i][j]=cM[i][j].imag();
        }
    }
    return M;
}

MapAnalyser::cfMatrix_t MapAnalyser::complexTypeMatrix_m(fMatrix_t M){
    cfMatrix_t cM;
    for (int i=0; i<6; i++){
        for (int j =0; j<6; j++){
            cM[i][j]=std::complex<double>(M[i][j],0);
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
            GSL_SET_COMPLEX(&temp,std::real(M[i][j]),std::imag(M[i][j]));
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

#if 0
//TODO Work in progress
void MapAnalyser::rearrangeEigen_m(cfMatrix_t& eigenVal, cfMatrix_t& EigenVec){
#ifdef PHIL_WRITE
    std::ofstream out;
    out.open("OUT.txt", std::ios::app);
#endif
    double precision = 1e-9;
    int pair = 0;
    std::complex<double> mult;
    std::pair<int, int> pairs[3];
    cfMatrix_t eigVa, eigVe, temp;

    out << "+++++++++++++++++++++++++++++++" << std::endl;
    //out << EigenVec << std::endl;
    for (int i = 0; pair < 3 && i <2*DIM; i++){
        for (int j = 0; pair < 3 && j < i; j++){
            mult=0;
            if (true){//std::abs(EigenVal[i][i].real() - EigenVal[j][j].real()) < precision
                    //&& std::abs(EigenVal[i][i].imag() + EigenVal[j][j].imag() ) < precision ){
                out << eigenVal[i][i] << "    " << eigenVal[j][j] << std::endl;
                out << "/////////////Calc///////////" << std::endl;
                for (int k = 0; k <DIM; k++){

                    out << -1. * EigenVec[2*k+1][i]  << EigenVec[2*k][j] << std::endl;
                    out << EigenVec[2*k][i] << EigenVec[2*k+1][j] << std::endl;


                    mult += -1. * EigenVec[2*k+1][i] * EigenVec[2*k][j];
                    mult += EigenVec[2*k][i] * EigenVec[2*k+1][j];
                    out <<"in between mult:  " <<  mult << std::endl;
                }

                out << mult << std::endl;
                out << "Criteria    " << std::abs(mult.imag()-1.) << std::endl;
                if (std::abs(mult.imag()-1.) < precision){
                    pairs[pair]= std::pair<int,int>(i,j);
                    out << i << "    " << j << std::endl;
                    pair++;
                }

            }

         }
    }
}
#endif
