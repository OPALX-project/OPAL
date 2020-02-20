#ifndef CLASSIC_PartPusher_H
#define CLASSIC_PartPusher_H

#include "Algorithms/Vektor.h"
#include "Algorithms/PartData.h"
#include "Physics/Physics.h"

/*


 */

class BorisPusher {
    
public:
    BorisPusher(const PartData &ref);
    BorisPusher();
    void initialise(const PartData *ref);
    
    void kick(const Vector_t &R, Vector_t &P,
              const Vector_t &Ef, const Vector_t &Bf,
              const double &dt) const;
    
    void kick(const Vector_t &R, Vector_t &P,
              const Vector_t &Ef, const Vector_t &Bf,
              const double &dt, const double &mass,
              const double &charge) const;
    
    void push(Vector_t &R, const Vector_t &P, const double &dt) const;

private:
    const PartData *itsReference;
};

inline BorisPusher::BorisPusher(const PartData &ref):
    itsReference(&ref)
{ }

inline BorisPusher::BorisPusher():
    itsReference(NULL)
{ }

inline void BorisPusher::initialise(const PartData *ref)
{
    itsReference = ref;
}

inline void BorisPusher::kick(const Vector_t &R, Vector_t &P,
                              const Vector_t &Ef, const Vector_t &Bf,
                              const double &dt) const
{
    kick(R, P, Ef, Bf, dt, itsReference->getM(), itsReference->getQ());
}


inline void BorisPusher::kick(const Vector_t &/*R*/, Vector_t &P,
                              const Vector_t &Ef, const Vector_t &Bf,
                              const double &dt, const double &mass,
                              const double &charge) const
{
    // Implementation follows chapter 4-4, p. 61 - 63 from
    // Birdsall, C. K. and Langdon, A. B. (1985). Plasma physics
    // via computer simulation.
    //
    // Up to finite precision effects, the new implementation is equivalent to the
    // old one, but uses less floating point operations.
    //
    // Relativistic variant implemented below is described in
    // chapter 15-4, p. 356 - 357.
    // However, since other units are used here, small
    // modifications are required. The relativistic variant can be derived
    // from the nonrelativistic one by replacing
    //     mass
    // by
    //     gamma * rest mass
    // and transforming the units.
    //
    // Parameters:
    //     R = x / (c * dt): Scaled position x, not used in here
    //     P = v / c * gamma: Scaled velocity v
    //     Ef: Electric field
    //     Bf: Magnetic field
    //     dt: Timestep
    //     mass = rest energy = rest mass * c * c
    //     charge

    using Physics::c;

    // Half step E
    P += 0.5 * dt * charge * c / mass * Ef;

    // Full step B

    /*
        LF
    double const gamma = sqrt(1.0 + dot(P, P));
    Vector_t const t = dt * charge * c * c / (gamma * mass) * Bf;
    P +=  cross(P, t);
    */

    double const gamma = sqrt(1.0 + dot(P, P));
    Vector_t const t = 0.5 * dt * charge * c * c / (gamma * mass) * Bf;
    Vector_t const w = P + cross(P, t);
    Vector_t const s = 2.0 / (1.0 + dot(t, t)) * t;
    P += cross(w, s);

    /* a poor Leap-Frog
        P += 1.0 * dt * charge * c * c / (gamma * mass) * cross(P,Bf);
    */

    // Half step E
    P += 0.5 * dt * charge * c / mass * Ef;
}


inline void BorisPusher::push(Vector_t &R, const Vector_t &P, const double &/* dt */) const {
    /** \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
     *
     * \code
     * R[i] += 0.5 * P[i] * recpgamma;
     * \endcode
     */
    R += 0.5 * P / sqrt(1.0 + dot(P, P));
}

#endif