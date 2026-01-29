//
// Class RK4
//   Fourth order Runge-Kutta time integrator
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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

template <typename FieldFunction, typename ... Arguments>
bool RK4<FieldFunction, Arguments ...>::doAdvance_m(PartBunchBase<double, 3>* bunch,
                                                    const size_t& i,
                                                    const double& t,
                                                    const double dt,
                                                    Arguments& ... args) const
{
    // Fourth order Runge-Kutta integrator
    // arguments:
    //   x          Current value of dependent variable
    //   t          Independent variable (usually time)
    //   dt         Step size (usually time step)
    //   i          index of particle
    
    double x[6];
    
    this->copyTo(bunch->R[i], bunch->P[i], &x[0]);
    
    double  deriv1[6];
    double  deriv2[6];
    double  deriv3[6];
    double  deriv4[6];
    double  xtemp[6];

    // Evaluate f1 = f(x,t).

    bool outOfBound = derivate_m(bunch, x, t, deriv1 , i, args ...);
    if (outOfBound)
        return false;

    // Evaluate f2 = f( x+dt*f1/2, t+dt/2 ).
    const double half_dt = 0.5 * dt;
    const double t_half = t + half_dt;

    for(int j = 0; j < 6; ++j)
        xtemp[j] = x[j] + half_dt * deriv1[j];

    outOfBound = derivate_m(bunch, xtemp, t_half, deriv2 , i, args ...);
    if (outOfBound)
        return false;

    // Evaluate f3 = f( x+dt*f2/2, t+dt/2 ).
    for(int j = 0; j < 6; ++j)
        xtemp[j] = x[j] + half_dt * deriv2[j];

    outOfBound = derivate_m(bunch, xtemp, t_half, deriv3 , i, args ...);
    if (outOfBound)
        return false;

    // Evaluate f4 = f( x+dt*f3, t+dt ).
    double t_full = t + dt;
    for(int j = 0; j < 6; ++j)
        xtemp[j] = x[j] + dt * deriv3[j];

    outOfBound = derivate_m(bunch, xtemp, t_full, deriv4 , i, args ...);
    if (outOfBound)
        return false;

    // Return x(t+dt) computed from fourth-order R-K.
    for(int j = 0; j < 6; ++j)
        x[j] += dt / 6.*(deriv1[j] + deriv4[j] + 2.*(deriv2[j] + deriv3[j]));
    
    this->copyFrom(bunch->R[i], bunch->P[i], &x[0]);

    return true;
}


template <typename FieldFunction, typename ... Arguments>
bool RK4<FieldFunction, Arguments ...>::derivate_m(PartBunchBase<double, 3>* bunch,
                                                   double* y,
                                                   const double& t,
                                                   double* yp,
                                                   const size_t& i,
                                                   Arguments& ... args) const
{
    // New for OPAL 2.0: Changing variables to m, T, s
    // Currently: m, ns, kG

    Vector_t externalE, externalB, tempR;

    externalB = Vector_t(0.0, 0.0, 0.0);
    externalE = Vector_t(0.0, 0.0, 0.0);

    for(int j = 0; j < 3; ++j)
        tempR(j) = y[j];

    bunch->R[i] = tempR;

    bool outOfBound = this->fieldfunc_m(t, i, externalE, externalB, args ...);

    double qtom = bunch->Q[i] / (bunch->M[i] * mass_coeff);   // m^2/s^2/GV

    double tempgamma = sqrt(1 + (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]));

    yp[0] = c_mtns / tempgamma * y[3]; // [m/ns]
    yp[1] = c_mtns / tempgamma * y[4]; // [m/ns]
    yp[2] = c_mtns / tempgamma * y[5]; // [m/ns]

/*
    yp[0] = c_mmtns / tempgamma * y[3]; // [mm/ns]
    yp[1] = c_mmtns / tempgamma * y[4]; // [mm/ns]
    yp[2] = c_mmtns / tempgamma * y[5]; // [mm/ns]
*/

    yp[3] = (externalE(0) / Physics::c  + (externalB(2) * y[4] - externalB(1) * y[5]) / tempgamma) * qtom; // [1/ns]
    yp[4] = (externalE(1) / Physics::c  - (externalB(2) * y[3] - externalB(0) * y[5]) / tempgamma) * qtom; // [1/ns];
    yp[5] = (externalE(2) / Physics::c  + (externalB(1) * y[3] - externalB(0) * y[4]) / tempgamma) * qtom; // [1/ns];

    return outOfBound;
}


template <typename FieldFunction, typename ... Arguments>
void RK4<FieldFunction, Arguments ...>::copyTo(const Vector_t& R,
                                               const Vector_t& P,
                                               double* x) const
{
    for (int j = 0; j < 3; j++) {
        x[j]   = R(j);    // [x,y,z] (mm)
        x[j+3] = P(j);  // [px,py,pz] (beta*gamma)
    }
}


template <typename FieldFunction, typename ... Arguments>
void RK4<FieldFunction, Arguments ...>::copyFrom(Vector_t& R,
                                                 Vector_t& P,
                                                 double* x) const
{
    for(int j = 0; j < 3; j++) {
        R(j) = x[j];    // [x,y,z] (mm)
        P(j) = x[j+3];  // [px,py,pz] (beta*gamma)
    }
}
