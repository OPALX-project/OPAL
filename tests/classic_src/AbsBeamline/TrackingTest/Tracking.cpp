#include <sstream>

#include "gsl/gsl_odeiv.h"
#include "gsl/gsl_errno.h"

#include "Utilities/GeneralClassicException.h"
#include "AbsBeamline/Component.h"
#include "Physics/Physics.h"

#include "classic_src/AbsBeamline/TrackingTest/Tracking.h"

const double Tracking::_c_l = 299.792458;
Tracking* Tracking::tracking = NULL;

void Tracking::integrate(Tracking::var eqm, double targetIndie,
                          double* y, Component* field, double step_size) {
    // y in ns, mm, MeV/c^2, MeV/c
    _m = ::sqrt(y[4]*y[4] - y[5]*y[5] - y[6]*y[6] - y[7]*y[7]);
    if( std::isnan(_m) || std::isinf(_m) ) {
        if (y[4]*y[4] - y[5]*y[5] - y[6]*y[6] - y[7]*y[7] > -1e-6) {
        _m = 0.;
        } else {
        throw(GeneralClassicException("Tracking::integrate", "Mass undefined in stepping function"));
        }
    }

    field_m = field;
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
    gsl_odeiv_step    * step    = gsl_odeiv_step_alloc(T,8);
    gsl_odeiv_control * control = NULL; //gsl_odeiv_control_y_new(_absoluteError, _relativeError);
    gsl_odeiv_evolve  * evolve  = gsl_odeiv_evolve_alloc(8);
    gsl_odeiv_system    system = {z_equations_motion, NULL, 8, NULL};
    if (eqm == Tracking::t) {
        system.function = t_equations_motion;
    }

    tracking = this;
    double h = step_size;
    int    nsteps = 0;
    double indie = 0.;
    while(fabs(indie-targetIndie) > 1e-6) {
        nsteps++;
        int status =  gsl_odeiv_evolve_apply(evolve, control, step, &system, &indie, targetIndie, &h, y);
        if(status != GSL_SUCCESS) {
            throw(GeneralClassicException("Tracking::integrate", "Failed during tracking"));
            break;
        }
        if(nsteps > _maxNSteps) {
            std::stringstream ios;
            ios << "Killing tracking with step size " << h << " at step " << nsteps << " of " << _maxNSteps << "\n" 
                << "t: " << y[0] << " pos: " << y[1] << " " << y[2] << " " << y[3] << "\n"
                << "E: " << y[4] << " mom: " << y[5] << " " << y[6] << " " << y[7] << std::endl; 
            throw(GeneralClassicException("Tracking::integrate", "Exceeded maximum number of stepsq"));
            break;
        }
    }
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
}


int Tracking::t_equations_motion(double /*t*/,   const double x[8], double dxdt[8], void* /*params*/)
{
    Vector_t R(x[1]/1e3, x[2]/1e3, x[3]/1e3), P, E, B; // mm
    tracking->field_m->apply(R, P, x[0], E, B);
    B *= 1e-4; // kT
    double q_c    = tracking->_q; // units of electron charge
    //dx/dt = px/E
    dxdt[1] = _c_l*x[5]/x[4];
    dxdt[2] = _c_l*x[6]/x[4];
    dxdt[3] = _c_l*x[7]/x[4];
    //dt/dt = 1.
    dxdt[0] = 1.;
    // dp/dt = qE + vxB
    dxdt[5] = q_c*_c_l*( dxdt[2]*B[2] - dxdt[3]*B[1] ) + q_c*E[0]*_c_l;
    dxdt[6] = q_c*_c_l*( dxdt[3]*B[0] - dxdt[1]*B[2] ) + q_c*E[1]*_c_l;
    dxdt[7] = q_c*_c_l*( dxdt[1]*B[1] - dxdt[2]*B[0] ) + q_c*E[2]*_c_l;
    // E^2  = px^2+py^2+pz^2+ m^2
    // E dE = px dpx + py dpy + pz dpz
    // dEdt = dxdt dpx + dydt dpy + dzdt dpz (ignore B as B conserves energy)
    dxdt[4] = q_c*E[0]*dxdt[1]+q_c*E[1]*dxdt[2]+q_c*E[2]*dxdt[3];
    return GSL_SUCCESS;
}

int Tracking::z_equations_motion(double /*z*/, const double x[8], double dxdz[8], void* /*params*/)
{
    if(::fabs(x[7]) < 1e-9) { return GSL_ERANGE;}
    Vector_t R(x[1], x[2], x[3]), P, E, B;
    tracking->field_m->apply(R, P, x[0], E, B);
    B *= 1e-4; // kT

    double q_c    = tracking->_q; // units of electron charge
    double dtdz   = x[4]/x[7]; //NOTE to self:- why not use x[4]/fabs(x[7]) to simplify things?
    double dir    = ::fabs(x[7])/x[7]; //direction of motion
    //dx/dz = px/pz
    dxdz[0] = dtdz/_c_l;
    dxdz[1] = x[5]/x[7];
    dxdz[2] = x[6]/x[7];
    dxdz[3] = 1.;
    // dp/dz = dp/dt dt/dz = dt/dz qE + dt/dz dx/dt X B
    dxdz[5] = q_c*_c_l*( dxdz[2]*B[2] - dxdz[3]*B[1] ) + q_c*E[0]*dtdz*dir; //dpx/dz = qc/pz py*bz
    dxdz[6] = q_c*_c_l*( dxdz[3]*B[0] - dxdz[1]*B[2] ) + q_c*E[1]*dtdz*dir;
    dxdz[7] = q_c*_c_l*( dxdz[1]*B[1] - dxdz[2]*B[0] ) + q_c*E[2]*dtdz*dir;
    // E^2  = px^2+py^2+pz^2+ m^2
    // E dE = px dpx + py dpy + pz dpz
    // dEdz = dxdt dpxdz + dydt dpydz + dzdt dpzdz (ignore B as B conserves energy)
    // dEdz = dxdz Ex + dydz Ey + Ez (ignore B as B conserves energy)
    dxdz[4] = (dxdz[1]*q_c*E[0]+dxdz[2]*q_c*E[1]+q_c*E[2])*dir;
    return GSL_SUCCESS;
}
