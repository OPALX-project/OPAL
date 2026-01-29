//
// Class: Hamiltonian
//   Constructs thick lens Hamiltonian up to arbitrary order for beamline elements
//
// Copyright (c) 2018, Philippe Ganz, ETH Zürich
// All rights reserved
//
// Implemented as part of the Master thesis
// "s-based maps from TPS & Lie-Series applied to Proton-Therapy Gantries"
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

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Classic/FixedAlgebra/FTps.h"
#include "Classic/FixedAlgebra/FVps.h"
#include "Classic/FixedAlgebra/FTpsMath.h"
#define PSdim 6

#include <functional>

class Hamiltonian
{

public:
    typedef FTps<double, PSdim> series_t;

    explicit Hamiltonian(int truncOrder);

    series_t x;     /*!< Phase space 1st dimension */
    series_t px;    /*!< Phase space 2nd dimension */
    series_t y;     /*!< Phase space 3rd dimension */
    series_t py;    /*!< Phase space 4th dimension */
    series_t z;     /*!< Phase space 5th dimension */
    series_t delta; /*!< Phase space 6th dimension */


    /**Drift Space Hamiltonian
     *
     * \f[
     * H_{Drift}= \frac{\delta}{\beta_0} -
     * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } }
     * \f]
     *
     * @param gamma0 Lorenz factor
     */
    Hamiltonian::series_t drift(const double& gamma0);

    /**:TODO: WIP:Rectangular Bend Hamiltonian
     * \f[
     * H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
     * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -p_x^2 -p_y^2 - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
     * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right)
     * \f]
    */
    Hamiltonian::series_t rbend(double& beta0,
                                double& gamma0,
                                double& q,
                                double& h,
                                double& K0);

    /**Sector Bend Hamiltonian
     *
     * \f[
     * H_{DipoleFringeField}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
     * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2
     * -\left( p_x - \frac{1}{2} \frac{k_0}{l} \left( s^2 - y^2 \right) \right) ^2 -p_y^2
     * - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } }
     * + \left( 1+ hx \right) \frac{1}{2} \frac{k_0}{l} \left( y^2 - x^2 \right) \tan \left( \Psi \right)
     * \f]
     *
     * @param gamma0 Lorenz factor
     * @param h curvature (\f$ \frac{1}{R} \f$ , where \f$ R \f$ is the bend radius)
     * @param k0 normalized magnetic field (\f$ k0 = \frac{B q}{P_0}\f$,
     * where \f$q\f$ is the particle charge and \f$P_0\f$ the momentum of the reference particle)
     */
    Hamiltonian::series_t sbend(const double& gamma0,
                                const double& h,
                                const double& k0);

    /**:TODO: WIP: Fringe Field SBend
     *
     * \f[H_{Dipole}= \frac{\delta}{\beta_0} - \left( 1+ hx \right)
     * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2 -(p_x - a_x)^2 -p_y^2
     * - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } } +
     * \left( 1+ hx \right) k_0 \left(x - \frac{hx^2}{2 \left( 1+ hx \right)}\right)
     * \f]
     *
     * @param beta0
     * @param gamma0 Lorenz factor
     * @param k0 normalized magnetic field (\f$ k0 = \frac{B q}{P_0}\f$,
     * where \f$q\f$ is the particle charge and \f$P_0\f$ the momentum of the reference particle)
     * @param ax Vector potential in x
     * @param az longitudinal vector potential (in z)
     */
    Hamiltonian::series_t bendFringe(double& beta0,
                                     double& gamma0,
                                     double& h,
                                     double& k0,
                                     series_t& ax,
                                     series_t& az);


    /**Quadrupole Hamiltonian
     *
     * \f[
     * H_{Quadrupole}= \frac{\delta}{\beta_0} -
     * \sqrt{\left(\frac{1}{\beta_0} + \delta \right)^2  -p_x^2 -p_y^2
     * - \frac{1}{\left(\beta_0 \gamma_0\right)^2 } }   +
     * \frac{1}{2} k_1 \left( x^2 - y^2 \right)
     * \f]
     *
     * @param gamma0 Lorenz factor
     * @param q particle charge
     * @param k1 normalised field gradient (\f$ k1 = \frac{B q}{P_0 r_0}\f$,
     * where \f$q\f$ is the particle charge, \f$P_0\f$ the momentum of the reference particle
     * and \f$r_0\f$ the element aperture)
     */
    Hamiltonian::series_t quadrupole(const double& gamma0,
                                     const double& q,
                                     const double& k1);


    /**Hamiltonian for a linear Thin Lens fringe field approximation
     * \f[
     * H_{ThinLens} = \frac{1}{2} (x^2 - y^2) k_0 \tan \left( \Psi \right)
     * \f]
     *
     * @param phi pole face roation angle
     * @param k0 normalized magnetic field (\f$ k0 = \frac{B q}{P_0}\f$,
     * where \f$q\f$ is the particle charge and \f$P_0\f$ the momentum of the reference particle)
     */
    Hamiltonian::series_t fringeField(const double& phi,
                                      const double& k0);

private:
    /**
     * Calculate \f[\beta = \sqrt{1-\frac{1}{\gamma^2}}
     *
     * @param gamma Lorenz factor
     */
    double getBeta_m(const double& gamma);

private:
    int truncOrder_m;  ///< Map truncation order (= Hamiltonian truncOrd - 1)

};


// /**Hamiltonian for Space Charges
//  * \f[ K = \frac{3 q I \lambda}{20 \sqrt{5} \pi \epsilon_0 m c^3 \beta^2 \gamma^3} , \;
//  *  H_{sc} = -\frac{1}{2} K_x  x^2
//  *  -\frac{1}{2} K_y  y^2
//  *  -\frac{1}{2} K_z  z^2 \gamma^2\f]
//  */
// void spaceCharge(series_t& H,
//                  PartBunchBase<double, 3>* bunch)
// {
//     //Derived form Distribution/SigmaGenerator.cpp
//     // convert m from MeV/c^2 to eV*s^{2}/m^{2}
//
//     double m = bunch->getM() / (Physics::c * Physics::c);
//
//     // formula (57)
//     double gamma = bunch->get_gamma();
//     double beta = std::sqrt( 1. - 1/( gamma*gamma ) );
//
//     double gamma2=gamma * gamma;
//
//     double freq = itsBeam_m->getFrequency() * 1e6;              // [freq] = Hz (form MHz)
//     double I = itsBeam_m->getCurrent();                         // [I] = A
//     double q = bunch->getQ();                              // [q] = e
//
// #ifdef PHIL_WRITE
//     std::ofstream dAngle;
//     dAngle.open("dipoleAngle.txt", std::ios::app);
// #endif
//
//     dAngle << "Freq:" <<  freq << std::endl;
//     dAngle << "I: " <<  I << std::endl;
//     dAngle << "SigmaMatirx" << bunch->getSigmaMatrix() << std::endl;
//     dAngle << "mass" << bunch->getM()<< std::endl;
//     dAngle << "charge" << q << std::endl;
//
//     //TODO check formula for lambda
//     double lam = Physics::c / freq;                              // wavelength, [lam] = m
//     double K3 = 3.0 * I * q * lam / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * m *
//                     Physics::c * Physics::c * Physics::c * beta * beta * gamma * gamma2);            // [K3] = m
//     dAngle << "K3:  " << K3 << std::endl;
//     dAngle << "K3:  " << lam << std::endl;
//     dAngle << "gamma:  " << gamma << std::endl;
//     dAngle << "beta:  " << beta << std::endl;
//     // formula (30), (31)
//     fMatrix_t sigmMatrix = bunch->getSigmaMatrix();
//
//     double sx = std::sqrt(std::fabs(sigmMatrix(0,0)));              //[sx] = m
//     double sy = std::sqrt(std::fabs(sigmMatrix(2,2)));              //[sy] = m
//     double sz = std::sqrt(std::fabs(sigmMatrix(4,4)));              //[sz] = m
//
//     dAngle << "sx:  " << sx << std::endl;
//     dAngle << "sz:  " << sz << std::endl;
//     double tmp = sx * sy;                                           // [tmp] = m^{2}
//
//     double f = std::sqrt(tmp) / (3.0 * gamma * sz);                 // [f] = 1
//     double kxy = K3 * std::fabs(1.0 - f) / ((sx + sy) * sz);        // [kxy] = 1/m
//
//     dAngle << "f:  " << f << std::endl;
//     dAngle << "kxy:  " << kxy << std::endl;
//     double Kx = kxy / sx;                                           //[Kx] = 1/m^{2}
//     double Ky = kxy / sy;                                           //[Ky] = 1/m^{2}
//     double Kz = K3 * f / (tmp * sz);                                //[Kz] = 1/m^{2}
//
//     // Change units for application on [m], [p/p0] and [E/p0-1/b0]
//     // (these units are for [mm], [mrad] and [promille])
//     // according Hinterberger - Physik der Teilcherbeschleuniger Eq.4.8
//
//     double mm2m = 1e3;
//     double mrad2pnorm= 1e3; // for small angles (tan(a) ~ a)
//     double promille2delta = 1/(10 * bunch->getInitialBeta()); // [E/p0-1/b0] = [1/b0 * (E0 - E)/E0]
//
//
//     H = ( -0.5 * Kx * x * x
//           -0.5 * Ky * y * y ) * mm2m / mrad2pnorm
//         - (0.5 * Kz * z * z * gamma2) * promille2delta ;
// }


#endif
