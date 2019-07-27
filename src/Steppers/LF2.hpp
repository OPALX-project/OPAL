#include "BorisPusher.h"

template <typename FieldFunction, typename ... Arguments>
bool LF2<FieldFunction, Arguments ...>::doAdvance_m(PartBunchBase<double, 3>* bunch,
                                                    const size_t& i,
                                                    const double& t,
                                                    const double dt,
                                                    Arguments& ... args) const
{
    bool flagNoDeletion = true;

    // push for first LF2 half step
    push_m(bunch->R[i], bunch->P[i], 0.5 * dt * 1.0e-9);  // ns --> s

    flagNoDeletion = kick_m(bunch, i, t, dt * 1.0e-9, args ...);

    // push for second LF2 half step
    push_m(bunch->R[i], bunch->P[i], 0.5 * dt * 1.0e-9);  // ns --> s

    return flagNoDeletion;
}


template <typename FieldFunction, typename ... Arguments>
void LF2<FieldFunction, Arguments ...>::push_m(Vector_t& R, const Vector_t& P,
                                               const double& h) const
{
    double const gamma = sqrt(1.0 + dot(P, P));
    double const c_gamma = Physics::c / gamma;
    Vector_t const v = P * c_gamma;
    R += h * v;
}


template <typename FieldFunction, typename ... Arguments>
bool LF2<FieldFunction, Arguments ...>::kick_m(PartBunchBase<double, 3>* bunch, const size_t& i,
                                               const double& t, const double& h,
                                               Arguments& ... args) const
{
    Vector_t externalE = Vector_t(0.0, 0.0, 0.0);
    Vector_t externalB = Vector_t(0.0, 0.0, 0.0);
    
    bool outOfBound = this->fieldfunc_m(t, i, externalE, externalB, args ...);
    
    if ( outOfBound )
        return false;
    
    
    double const q = bunch->Q[0] / Physics::q_e; // For now all particles have the same charge
    double const M = bunch->M[0] * 1.0e9; // For now all particles have the same rest energy
    
    BorisPusher pusher;
    
    pusher.kick(bunch->R[i], bunch->P[i],
                externalE, externalB,
                h, M, q);
    
    return true;
}
