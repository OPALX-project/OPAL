//
// Class Harmonics
//   This class computes the cyclotron map based on harmonics.
//   All functions are copied and translated to C++ from the original program inj2_ana.c of Dr. C. Baumgarten.
//
// Copyright (c) 2014 - 2015, Christian Baumgarten, Paul Scherrer Institut, Villigen PSI, Switzerland
//                            Matthias Frey, ETH Zürich
// All rights reserved
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
#ifndef HARMONICS_H
#define HARMONICS_H

#include <cmath>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "Physics/Physics.h"

#include <boost/numeric/ublas/matrix.hpp>

#include "matrix_vector_operation.h"

template<typename Value_type, typename Size_type>
class Harmonics
{
public:
    /// Type of variable
    typedef Value_type value_type;
    /// Type of size
    typedef Size_type size_type;
    /// Type for specifying matrices
    typedef boost::numeric::ublas::matrix<value_type> matrix_type;


    /// Defines starting point of computation
    enum START { SECTOR_ENTRANCE, SECTOR_CENTER, SECTOR_EXIT, VALLEY_CENTER };

    /// Constructor
    /*!
     * @param wo is the nominal orbital frequency in Hz
     * @param Emin is the minimal energy in MeV
     * @param Emax is the maximal energy in MeV
     * @param nth is the number of angle steps
     * @param nr is the number of radial steps
     * @param nSector is the number of sectors
     * @param E is the kinetic energy [MeV]
     * @param E0 is the potential energy [MeV]
     */
    Harmonics(value_type, value_type, value_type, size_type, size_type, size_type, value_type, value_type);

    /// Compute all maps of the cyclotron using harmonics
    /*!
     * @param filename1 name of file 1
     * @param filename2 name of file 2
     * @param startmod (in center of valley, start of sector, ...)
     */
    std::vector<matrix_type> computeMap(const std::string, const std::string, const int);

    /// Returns the radial and vertical tunes
    std::pair<value_type, value_type> getTunes();

    /// Returns the radius
    value_type getRadius();

    /// Returns step size
    value_type getPathLength();

private:

    /// Stores the nominal orbital frequency in Hz
    value_type wo_m;
    /// Stores the minimum energy in MeV
    value_type Emin_m;
    /// Stores the maximum energy in MeV
    value_type Emax_m;
    /// Stores the number of angle splits
    size_type nth_m;
    /// Stores the number of radial splits
    size_type nr_m;
    /// Stores the number of sectors
    size_type nSector_m;
    /// Stores the tunes (radial, vertical)
    std::pair<value_type, value_type> tunes_m;
    /// Stores the radius
    value_type radius_m;
    /// Stores the path length
    value_type ds_m;
    /// Stores the energy for which we perform the computation
    value_type E_m;
    /// Stores the potential energy [MeV]
    value_type E0_m;

    /// Compute some matrix (ask Dr. C. Baumgarten)
    matrix_type __Mix6(value_type, value_type, value_type);
    /// Compute some matrix  (ask Dr. C. Baumgarten)
    matrix_type __Md6(value_type, value_type);
    /// Compute some matrix  (ask Dr. C. Baumgarten)
    matrix_type __Mb6(value_type, value_type, value_type);
    /// Compute some matrix  (ask Dr. C. Baumgarten)
    matrix_type __Mb6k(value_type, value_type, value_type, value_type);
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
Harmonics<Value_type, Size_type>::Harmonics(value_type wo, value_type Emin, value_type Emax,
    size_type nth, size_type nr, size_type nSector, value_type E, value_type E0)
    : wo_m(wo), Emin_m(Emin), Emax_m(Emax), nth_m(nth), nr_m(nr), nSector_m(nSector), E_m(E), E0_m(E0)
{}

template<typename Value_type, typename Size_type>
std::vector<typename Harmonics<Value_type, Size_type>::matrix_type> Harmonics<Value_type, Size_type>::computeMap(const std::string filename1, const std::string filename2, const int startmod) {
    // i --> different energies
    // j --> different angles

    // ---------------


//     unsigned int nth_m = 360; //270;
    unsigned int N   = 4;
    value_type two_pi    = 2.0 * M_PI;
    value_type tpiN      = two_pi / value_type(N);
    value_type piN       = M_PI / value_type(N);
//     unsigned int nr_m  = 180;
    unsigned int cnt = 0;
    bool set = false;

    value_type gap = 0.05;
    value_type K1  = 0.45;

    std::vector<matrix_type> Mcyc(nth_m);

    value_type beta_m;

    std::vector<value_type> gamma(nr_m), E(nr_m), PC(nr_m), nur(nr_m), nuz(nr_m);

    std::vector<value_type> R(nr_m), r(nr_m);
    std::vector<value_type> Bmag(nr_m), k(nr_m);
    std::vector<value_type> alpha(nr_m), beta(nr_m);

    //----------------------------------------------------------------------
    // read files
    std::ifstream infile2(filename2);
    std::ifstream infile1(filename1);

    if (!infile1.is_open() || !infile2.is_open()) {
        std::cerr << "Couldn't open files!" << std::endl;
        std::exit(0);
    }

    unsigned int n = 0;
    value_type rN1, cN1, sN1, aN1, phN1, rN2, cN2, sN2, aN2, phN2;

    while (infile1 >> rN1 >> cN1 >> sN1 >> aN1 >> phN1 &&
           infile2 >> rN2 >> cN2 >> sN2 >> aN2 >> phN2) {

        R[n]      = rN1 * 0.01;
        alpha[n]  = 2.0 * std::acos(aN2 / aN1) / value_type(N);
        beta[n]   = std::atan(sN1 / cN1) / value_type(N);
        value_type f4 = aN1 * M_PI * 0.5 / std::sin(0.5 * alpha[n] * value_type(N));
        value_type f8 = aN2 * M_PI / std::sin(alpha[n] * value_type(N));

        Bmag[n] = 0.5 * (f4 + f8);

        n++;
    }

    infile1.close();
    infile2.close();
    //----------------------------------------------------------------------

    value_type s = std::sin(piN);
    value_type c = std::cos(piN);

    std::vector<value_type> dalp(nr_m), dbet(nr_m), len2(nr_m), len(nr_m), t1eff(nr_m), t2eff(nr_m);
    std::vector<value_type> t1(nr_m), t2(nr_m), t3(nr_m);

    dalp[0]   = (alpha[1]   - alpha[0])   / (R[1] - R[0]);
    dalp[n-1] = (alpha[n-1] - alpha[n-2]) / (R[n-1] - R[n-2]);
    dbet[0]   = (beta[1]    - alpha[0])   / (R[1] - R[0]);
    dbet[n-1] = (beta[n-1]  - beta[n-2])  / (R[n-1] - R[n-2]);

    for (unsigned int i = 1; i < n - 1; ++i) {
        dalp[i] = R[i] * (alpha[i+1] - alpha[i-1]) / (R[i+1] - R[i-1]);
        dbet[i] = R[i] * (beta[i+1]  - beta[i-1])  / (R[i+1] - R[i-1]);
    }

    for (unsigned int i = 0; i < n; ++i) {
        value_type cot = 1.0 / std::tan(0.5 * alpha[i]);

        value_type eps1 = std::atan(dbet[i] - 0.5 * dalp[i]);
        value_type eps2 = std::atan(dbet[i] + 0.5 * dalp[i]);

        value_type phi = piN - 0.5 * alpha[i];

        // bending radius
        r[i]    = R[i] * sin(0.5 * alpha[i]) / s;
        len2[i] = R[i] * std::sin(phi);
        len[i]  = two_pi * r[i] + 2.0 * len2[i] * value_type(N);

        value_type g1 =   phi + eps1;
        value_type g2 = - phi + eps2;

        t1[i] = std::tan(g1);
        t2[i] = std::tan(g2);
        t3[i] = std::tan(phi);

        value_type fac = gap / r[i] * K1;
        value_type psi = fac * (1.0 + std::sin(g1) * std::sin(g1)) / std::cos(g1);

        t1eff[i] = std::tan(g1 - psi);

        psi = fac * (1.0 + std::sin(g2) * std::sin(g2)) / std::cos(g2);
        t2eff[i] = std::tan(g2 + psi);

        beta_m     = wo_m / Physics::c * len[i] / two_pi;
        gamma[i] = 1.0 / std::sqrt(1.0 - beta_m * beta_m);
        E[i]     = E0_m * (gamma[i] - 1.0);
        PC[i]    = E0_m * gamma[i] * beta_m;
        Bmag[i]  = E0_m * 1.0e6 / Physics::c * beta_m * gamma[i] / r[i] * 10.0;

        if (!set && E[i] >= E_m) {
            set = true;
            cnt = i;
        }
    }

    // (average) field gradient
    for (unsigned int i = 0; i < n; ++i) {
        value_type dBdr;
        if (i == 0)
            dBdr = (Bmag[1] - Bmag[0]) / (R[1] - R[0]);
        else if (i == n - 1)
            dBdr = (Bmag[n-1] - Bmag[n-2]) / (R[n-1] - R[n-2]);
        else
            dBdr = (Bmag[i+1] - Bmag[i-1]) / (R[i+1] - R[i-1]);

        value_type bb  = std::sin(piN - 0.5 * alpha[i]) / std::sin(0.5 * alpha[i]);
        value_type eps = 2.0 * bb / (1.0 + bb * bb);
        value_type BaV = Bmag[i];

        k[i] = r[i] * r[i] / (R[i] * BaV) * dBdr * (1.0 + bb * value_type(N) / M_PI * s);
    }

    for (unsigned int i = 0; i < n; ++i) {
        if ((k[1] > -0.999) && (E[i] > Emin_m) && (E[i] < Emax_m + 1.0)) {
            value_type kx = std::sqrt(1.0 + k[i]);
            value_type Cx = std::cos(tpiN * kx);
            value_type Sx = std::sin(tpiN * kx);
            value_type pm, ky, Cy, Sy;

            if (k[i] >= 0.0) {
                pm = 1.0;
                ky = std::sqrt(std::abs(k[i]));
                Cy = std::cosh(tpiN * ky);
                Sy = std::sinh(tpiN * ky);
            } else {
                pm = -1.0;
                ky = std::sqrt(std::abs(k[i]));
                Cy = std::cos(tpiN * ky);
                Sy = std::sin(tpiN * ky);
            }

            value_type tav  = 0.5 * (t1[i] - t2[i]);
            value_type tmul = t1[i] * t2[i];
            value_type lrho = 2.0 * len2[i] / r[i];

            nur[i] = std::acos(Cx * (1.0 + lrho * tav) + Sx / kx * (tav - lrho * 0.5 * (kx * kx + tmul))) / tpiN;

            tav  = 0.5 * (t1eff[i] - t2eff[i]);
            tmul = t1eff[i] * t2eff[i];

            nuz[i] = std::acos(Cy * (1.0 - lrho * tav) + Sy / ky * (0.5 * lrho * (pm * ky * ky - tmul) - tav)) / tpiN;

//             std::cout << E[i] << " " << R[i] << " " << nur[i] << " " << nuz[i] << std::endl;
        }
    }
    // -------------------------------------------------------------------------

    int i = cnt;

    tunes_m = { nur[i], nuz[i] };
    radius_m = R[i];

//     for (int i = 9; i < 10; ++i) {       // n = 1 --> only for one energy
        value_type smax, slim, ss, gam2;

        smax = len[i] / value_type(N);
        slim = r[i] * tpiN;
        ds_m   = smax / value_type(nth_m);
        ss = 0.0;
        gam2 = gamma[i] * gamma[i];

        matrix_type Mi = __Mix6(  t1[i],   t1eff[i], r[i]);
        matrix_type Mx = __Mix6(- t2[i], - t2eff[i], r[i]);
        matrix_type Md = __Md6(ds_m, gam2);
        matrix_type Mb = __Mb6k(ds_m / r[i], r[i], k[i], gam2);

        unsigned int j = 0;

        matrix_type M0 = boost::numeric::ublas::zero_matrix<value_type>(6);
        matrix_type M1 = boost::numeric::ublas::zero_matrix<value_type>(6);

        switch (startmod) {
            case 0: // SECTOR_ENTRANCE
                ss = slim - ds_m;
                Mcyc[j++] = boost::numeric::ublas::prod(Mb, Mi);

                while (ss > ds_m) {
                    Mcyc[j++] = Mb;
                    ss -= ds_m;
                }

                M0 = __Mb6k(ss / r[i], r[i], k[i], gam2);
                M1 = __Md6(ds_m - ss, gam2);

                Mcyc[j++] = matt_boost::gemmm<matrix_type>(M1, Mx, M0);

                ss = smax - slim - (ds_m - ss);

                while ((ss > ds_m) && (j < nth_m)) {
                    Mcyc[j++] = Md;
                    ss -= ds_m;
                }

                if (j < nth_m)
                    Mcyc[j++] = __Md6(ss, gam2);

                break;

            case 1: // SECTOR_CENTER
                j  = 0;
                ss = slim * 0.5;

                while (ss > ds_m) {
                    Mcyc[j++] = Mb;
                    ss -= ds_m;
                }

                M0 = __Mb6k(ss / r[i], r[i], k[i], gam2);
                M1 = __Md6(ds_m - ss, gam2);

                Mcyc[j++] = matt_boost::gemmm<matrix_type>(M1, Mx, M0);

                ss = smax - slim - (ds_m - ss);

                while (ss > ds_m) {
                    Mcyc[j++] = Md;
                    ss -= ds_m;
                }

                M0 = __Md6(ss, gam2);
                M1 = __Mb6k((ds_m - ss) / r[i], r[i], k[i], gam2);

                Mcyc[j++] = matt_boost::gemmm<matrix_type>(M1, Mi, M0);

                ss = slim * 0.5 - (ds_m - ss);

                while ((ss > ds_m) && (j < nth_m)) {
                    Mcyc[j++] = Mb;
                    ss -= ds_m;
                }

                if (j < nth_m)
                    Mcyc[j++] = __Mb6k(ss / r[i], r[i], k[i], gam2);

                break;

            case 2: // SECTOR_EXIT
                j = 0;

                Mcyc[j++] = boost::numeric::ublas::prod(Md,Mx);

                ss = smax - slim - ds_m;

                while (ss > ds_m) {
                    Mcyc[j++] = Md;
                    ss -= ds_m;
                }

                M0 = __Md6(ss, gam2);
                M1 = __Mb6k((ds_m - ss) / r[i], r[i], k[i], gam2);

                Mcyc[j++] = matt_boost::gemmm<matrix_type>(M1, Mi, M0);

                ss = slim - (ds_m - ss);

                while ((ss > ds_m) && (j < nth_m)) {
                    Mcyc[j++] = Mb;
                    ss -= ds_m;
                }

                if (j < nth_m)
                    Mcyc[j++] = __Mb6k(ss / r[i], r[i], k[i], gam2);

                break;

            case 3: // VALLEY_CENTER
            default:
                j = 0;
                ss = (smax - slim) * 0.5;

                while (ss > ds_m) {
                    Mcyc[j++] = Md;
                    ss -= ds_m;
                }

                M0 = __Md6(ss, gam2);
                M1 = __Mb6k((ds_m - ss) / r[i], r[i], k[i], gam2);

                Mcyc[j++] = matt_boost::gemmm<matrix_type>(M1, Mi, M0);

                ss = slim - (ds_m - ss);
                while (ss > ds_m) {
                    Mcyc[j++] = Mb;
                    ss -= ds_m;
                }

                M0 = __Mb6k(ss / r[i], r[i], k[i], gam2);
                M1 = __Md6(ds_m - ss, gam2);

                Mcyc[j++] = matt_boost::gemmm<matrix_type>(M1, Mx, M0);

                ss = (smax - slim) * 0.5 - (ds_m - ss);
                while ((ss > ds_m) && (j < nth_m)) {
                    Mcyc[j++] = Md;
                    ss -= ds_m;
                }

                if (j < nth_m)
                    Mcyc[j++] = __Md6(ss, gam2);

                break;
        } // end of switch
//     } // end of for

    return std::vector<matrix_type>(Mcyc);
}

template<typename Value_type, typename Size_type>
inline std::pair<Value_type, Value_type> Harmonics<Value_type, Size_type>::getTunes() {
    return tunes_m;
}

template<typename Value_type, typename Size_type>
inline typename Harmonics<Value_type, Size_type>::value_type Harmonics<Value_type, Size_type>::getRadius() {
    return radius_m;
}

template<typename Value_type, typename Size_type>
inline typename Harmonics<Value_type, Size_type>::value_type Harmonics<Value_type, Size_type>::getPathLength() {
    return ds_m;
}


// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
typename Harmonics<Value_type, Size_type>::matrix_type Harmonics<Value_type, Size_type>::__Mix6(value_type tx, value_type ty, value_type r) {
    matrix_type M = boost::numeric::ublas::identity_matrix<value_type>(6);

    M(1,0) =   tx / r;
    M(3,2) = - ty / r;

    return M;
}

template<typename Value_type, typename Size_type>
typename Harmonics<Value_type, Size_type>::matrix_type Harmonics<Value_type, Size_type>::__Md6(value_type l, value_type gam2) {
    matrix_type M = boost::numeric::ublas::identity_matrix<value_type>(6);

    M(0,1) = l;
    M(2,3) = l;
    M(4,5) = l / gam2;

    return M;
}

template<typename Value_type, typename Size_type>
typename Harmonics<Value_type, Size_type>::matrix_type Harmonics<Value_type, Size_type>::__Mb6(value_type phi, value_type r, value_type gam2) {
    matrix_type M = boost::numeric::ublas::identity_matrix<value_type>(6);

    value_type C = std::cos(phi);
    value_type S = std::sin(phi);
    value_type l =r * phi;

    M(0,0) = M(1,1) = C;
    M(0,1) =   S * r;
    M(1,0) = - S / r;
    M(2,3) = l;
    M(4,5) = l / gam2 - l + r * S;
    M(4,0) = - (M(1,5) = S);
    M(4,1) = - (M(0,5) = r * (1.0 - C));

    return M;
}

template<typename Value_type, typename Size_type>
typename Harmonics<Value_type, Size_type>::matrix_type Harmonics<Value_type, Size_type>::__Mb6k(value_type phi, value_type r, value_type k, value_type gam2) {
    matrix_type M = boost::numeric::ublas::identity_matrix<value_type>(6);

    value_type C,S;

    if (k == 0.0) return __Mb6(phi,r,gam2);

    value_type fx = std::sqrt(1.0 + k);
    value_type fy = std::sqrt(std::abs(k));
    value_type c  = std::cos(phi * fx);
    value_type s  = std::sin(phi * fx);
    value_type l  = r * phi;


    if (k > 0.0) {
        C = std::cosh(phi * fy);
        S = std::sinh(phi * fy);
    } else {
        C = std::cos(phi * fy);
        S = std::sin(phi * fy);
    }

    M(0,0) = M(1,1) = c;
    M(0,1) =   s * r / fx;
    M(1,0) = - s * fx / r;
    M(2,2) = M(3,3) = C;
    M(2,3) = S * r / fy;

    value_type sign = (std::signbit(k)) ? value_type(-1) : value_type(1);

    M(3,2) = sign * S * fy / r;
    M(4,5) = l / gam2 - r / (1.0 + k) * (phi - s / fx);
    M(4,0)= - (M(1,5) = s / fx);
    M(4,1)= - (M(0,5) = r * (1.0 - c) / (1.0 + k));

    return M;
}

#endif