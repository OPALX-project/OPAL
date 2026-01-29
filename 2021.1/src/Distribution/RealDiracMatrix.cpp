//
// Class: RealDiracMatrix
//   Real Dirac matrix class.
//   They're ordered after the paper of Dr. C. Baumgarten: "Use of real Dirac matrices in two-dimensional coupled linear optics".
//   The diagonalizing method is based on the paper "Geometrical method of decoupling" (2012) of Dr. C. Baumgarten.
//
// Copyright (c) 2014, 2020 Christian Baumgarten, Paul Scherrer Institut, Villigen PSI, Switzerland
//                          Matthias Frey, ETH Zürich
// All rights reserved
//
// Implemented as part of the Semester thesis by Matthias Frey
// "Matched Distributions in Cyclotrons"
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
#include "RealDiracMatrix.h"

#include "Utilities/OpalException.h"

#include <cmath>
#include <string>
#include "matrix_vector_operation.h"

RealDiracMatrix::RealDiracMatrix() {};

typename RealDiracMatrix::sparse_matrix_t
RealDiracMatrix::getRDM(short i) {
    sparse_matrix_t rdm(4,4,4); // #nrows, #ncols, #non-zeros
    switch(i) {
        case 0:  rdm(0,1) = rdm(2,3) = 1; rdm(1,0) = rdm(3,2) = -1; break;
        case 1:  rdm(0,1) = rdm(1,0) = -1; rdm(2,3) = rdm(3,2) = 1; break;
        case 2:  rdm(0,3) = rdm(1,2) = rdm(2,1) = rdm(3,0) = 1;     break;
        case 3:  rdm(0,0) = rdm(2,2) = -1; rdm(1,1) = rdm(3,3) = 1; break;
        case 4:  rdm(0,0) = rdm(3,3) = -1; rdm(1,1) = rdm(2,2) = 1; break;
        case 5:  rdm(0,2) = rdm(2,0) = 1; rdm(1,3) = rdm(3,1) = -1; break;
        case 6:  rdm(0,1) = rdm(1,0) = rdm(2,3) = rdm(3,2) = 1;     break;
        case 7:  rdm(0,3) = rdm(2,1) = 1; rdm(1,2) = rdm(3,0) = -1; break;
        case 8:  rdm(0,1) = rdm(3,2) = 1; rdm(1,0) = rdm(2,3) = -1; break;
        case 9:  rdm(0,2) = rdm(1,3) = -1; rdm(2,0) = rdm(3,1) = 1; break;
        case 10: rdm(0,2) = rdm(3,1) = 1; rdm(1,3) = rdm(2,0) = -1; break;
        case 11: rdm(0,2) = rdm(1,3) = rdm(2,0) = rdm(3,1) = -1;    break;
        case 12: rdm(0,0) = rdm(1,1) = -1; rdm(2,2) = rdm(3,3) = 1; break;
        case 13: rdm(0,3) = rdm(3,0) = -1; rdm(1,2) = rdm(2,1) = 1; break;
        case 14: rdm(0,3) = rdm(1,2) = -1; rdm(2,1) = rdm(3,0) = 1; break;
        case 15: rdm(0,0) = rdm(1,1) = rdm(2,2) = rdm(3,3) = 1;     break;
        default: throw OpalException("RealDiracMatrix::getRDM()",
                                     "Index (i = " + std::to_string(i)
                                     + " out of range: 0 <= i <= 15"); break;
    }
    return rdm;
}


void RealDiracMatrix::diagonalize(matrix_t& Ms, sparse_matrix_t& R, sparse_matrix_t& invR) {

    // R and invR store the total transformation
    R = boost::numeric::ublas::identity_matrix<double>(4);
    invR = R;

    vector_t P(3), E(3), B(3), b;
    double mr, mg, mb, eps;

    // Lambda function to compute vectors E, P, B and scalar eps (it takes the current Ms as reference argument (--> [&])
    auto mult = [&](short i) {
        /*
         * For computing E, P, B, eps according to formula (C4) from paper:
         * Geometrical method of decoupling
         */
        matrix_t tmp = boost::numeric::ublas::prod(Ms,getRDM(i))+boost::numeric::ublas::prod(getRDM(i),Ms);
        return 0.125*matt_boost::trace(tmp);
    };

    // 1. Transformation with \gamma_{0}
    P = vector_t({ mult(1),  mult(2),  mult(3)});
    E = vector_t({ mult(4),  mult(5),  mult(6)});
    B = vector_t({-mult(7), -mult(8), -mult(9)});
    mr = boost::numeric::ublas::inner_prod(E, B);        // formula (31), paper: Geometrical method of decoupling
    mg = boost::numeric::ublas::inner_prod(B, P);        // formula (31), paper: Geometrical method of decoupling

    transform(Ms, 0, 0.5 * std::atan2(mg,mr), R, invR);

    // 2. Transformation with \gamma_{7}
    eps =  -mult(0);
    P   = vector_t({ mult(1),  mult(2),  mult(3)});
    E   = vector_t({ mult(4),  mult(5),  mult(6)});
    B   = vector_t({-mult(7), -mult(8), -mult(9)});
    b   = eps * B + matt_boost::cross_prod(E, P);      // formula (32), paper: Geometrical method of decoupling

    transform(Ms, 7, 0.5 * std::atan2(b(2), b(1)), R, invR);

    // 3. Transformation with \gamma_{9}
    eps =  -mult(0);
    P   = vector_t({ mult(1),  mult(2),  mult(3)});
    E   = vector_t({ mult(4),  mult(5),  mult(6)});
    B   = vector_t({-mult(7), -mult(8), -mult(9)});
    b   = eps * B + matt_boost::cross_prod(E, P);

    transform(Ms, 9, -0.5 * std::atan2(b(0), b(1)), R, invR);

    // 4. Transformation with \gamma_{2}
    eps =  -mult(0);
    P   = vector_t({ mult(1),  mult(2),  mult(3)});
    E   = vector_t({ mult(4),  mult(5),  mult(6)});
    B   = vector_t({-mult(7), -mult(8), -mult(9)});
    mr  = boost::numeric::ublas::inner_prod(E, B);
    b   = eps * B + matt_boost::cross_prod(E, P);

    // Transformation distinction made according to function "rdm_Decouple_F"
    // in rdm.c of Dr. Christian Baumgarten

    if (std::fabs(mr) < std::fabs(b(1))) {
        transform(Ms, 2, 0.5 * std::atanh(mr / b(1)), R, invR);
    } else {
        transform(Ms, 2, 0.5 * std::atanh(b(1) / mr), R, invR);
    }

    eps =  -mult(0);
    P   = vector_t({ mult(1),  mult(2),  mult(3)});
    E   = vector_t({ mult(4),  mult(5),  mult(6)});
    B   = vector_t({-mult(7), -mult(8), -mult(9)});

    // formula (31), paper: Geometrical method of decoupling
    mr = boost::numeric::ublas::inner_prod(E, B);
    mg = boost::numeric::ublas::inner_prod(B, P);
    mb = boost::numeric::ublas::inner_prod(E, P);

    double P2 = boost::numeric::ublas::inner_prod(P, P);
    double E2 = boost::numeric::ublas::inner_prod(E, E);

    // 5. Transform with \gamma_{0}
    transform(Ms, 0, 0.25 * std::atan2(mb,0.5 * (E2 - P2)), R, invR);

    // 6. Transformation with \gamma_{8}
    P(0) = mult(1); P(2) = mult(3);

    transform(Ms, 8, -0.5 * std::atan2(P(2), P(0)), R, invR);
}


RealDiracMatrix::sparse_matrix_t
RealDiracMatrix::diagonalize(matrix_t& sigma) {
    matrix_t S = boost::numeric::ublas::prod(sigma, getRDM(0));

    sparse_matrix_t R     = boost::numeric::ublas::identity_matrix<double>(4);
    sparse_matrix_t invR  = boost::numeric::ublas::identity_matrix<double>(4);
    sparse_matrix_t iRtot = boost::numeric::ublas::identity_matrix<double>(4);

    for (int i = 0; i < 6; ++i) {
        update(sigma, i, R, invR);

        iRtot = boost::numeric::ublas::prod(iRtot, invR);
        S = matt_boost::gemmm<matrix_t>(R, S, invR);
        sigma = - boost::numeric::ublas::prod(S, getRDM(0));
    }

    return iRtot;
}


typename RealDiracMatrix::matrix_t
RealDiracMatrix::symplex(const matrix_t& M) {
    /*
     * formula(16), p. 3
     */
    sparse_matrix_t rdm0 = getRDM(0);
    return 0.5 * (M + matt_boost::gemmm<matrix_t>(rdm0,boost::numeric::ublas::trans(M),rdm0));
}


void RealDiracMatrix::transform(matrix_t& M, short i, double phi,
                                sparse_matrix_t& Rtot, sparse_matrix_t& invRtot)
{
    if (phi) {  // if phi == 0 --> nothing happens, since R and invR would be identity_matrix matrix
        sparse_matrix_t R(4,4), invR(4,4);

        transform(i, phi, R, invR);

        // update matrices
        M = matt_boost::gemmm<matrix_t>(R,M,invR);
        Rtot = boost::numeric::ublas::prod(R,Rtot);
        invRtot = boost::numeric::ublas::prod(invRtot,invR);
    }
}


void RealDiracMatrix::transform(short i, double phi,
                                sparse_matrix_t& R, sparse_matrix_t& invR)
{
    if (phi) {  // if phi == 0 --> nothing happens, since R and invR would be identity_matrix matrix
        sparse_matrix_t I = boost::numeric::ublas::identity_matrix<double>(4);

        if ((i < 7 && i != 0) || (i > 10 && i != 14)) {
            R = I * std::cosh(phi) + getRDM(i) * std::sinh(phi);
            invR = I * std::cosh(phi) - getRDM(i) * std::sinh(phi);
        } else {
            R = I * std::cos(phi) + getRDM(i) * std::sin(phi);
            invR = I * std::cos(phi) - getRDM(i) * std::sin(phi);
        }
    }
}


void RealDiracMatrix::update(matrix_t& sigma, short i, sparse_matrix_t& R,
                             sparse_matrix_t& invR)
{
    double s0 =  0.25 * ( sigma(0, 0) + sigma(1, 1) + sigma(2, 2) + sigma(3, 3) );
    double s1 =  0.25 * (-sigma(0, 0) + sigma(1, 1) + sigma(2, 2) - sigma(3, 3) );
    double s2 =  0.5  * ( sigma(0, 2) - sigma(1, 3) );
    double s3 =  0.5  * ( sigma(0, 1) + sigma(2, 3) );
    double s4 =  0.5  * ( sigma(0, 1) - sigma(2, 3) );
    double s5 = -0.5  * ( sigma(0, 3) + sigma(1, 2) );
    double s6 =  0.25 * ( sigma(0, 0) - sigma(1, 1) + sigma(2, 2) - sigma(3, 3) );
    double s7 =  0.5  * ( sigma(0, 2) + sigma(1, 3) );
    double s8 =  0.25 * ( sigma(0, 0) + sigma(1, 1) - sigma(2, 2) - sigma(3, 3) );
    double s9 =  0.5  * ( sigma(0, 3) - sigma(1, 2) );

    vector_t P(3); P(0) = s1; P(1) = s2; P(2) = s3;
    vector_t E(3); E(0) = s4; E(1) = s5; E(2) = s6;
    vector_t B(3); B(0) = s7; B(1) = s8; B(2) = s9;

    double mr = boost::numeric::ublas::inner_prod(E, B);
    double mg = boost::numeric::ublas::inner_prod(B, P);
    double mb = boost::numeric::ublas::inner_prod(E, P);

    vector_t b = -s0 * B + matt_boost::cross_prod(E, P);

    switch (i) {
        case 0:
        {
            double eps = std::atan2(mg, mr);
            transform(0, 0.5 * eps, R, invR);
            break;
        }
        case 1:
        {
            double eps = std::atan2(b(2), b(1));
            transform(7, 0.5 * eps, R, invR);
            break;
        }
        case 2:
        {
            double eps = - std::atan2(b(0), b(1));
            transform(9, 0.5 * eps, R, invR);
            break;
        }
        case 3:
        {
            double eps = - std::atanh(mr / b(1));
            transform(2, 0.5 * eps, R, invR);
            break;
        }
        case 4:
        {
            double eps = 0.5 * std::atan2(2.0 * mb,
                                          boost::numeric::ublas::inner_prod(E, E) -
                                          boost::numeric::ublas::inner_prod(P, P));
            transform(0, 0.5 * eps, R, invR);
            break;
        }
        case 5:
        {
            double eps = - std::atan2(P(2), P(0));
            transform(8, 0.5 * eps, R, invR);
            break;
        }
        default:
        {
            throw OpalException("RealDiracMatrix::update()",
                                "Case " + std::to_string(i) +
                                " not available.");
        }
    }
}