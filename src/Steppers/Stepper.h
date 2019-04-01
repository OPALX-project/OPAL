#ifndef STEPPER_H
#define STEPPER_H

#include "Algorithms/PartBunchBase.h"
#include "Algorithms/Vektor.h"

#include <functional>

/*!
 * @precondition The field function has to return a
 * boolean and take at least the following arguments in that order:
 *  - double    specifying the time
 *  - int       specifying the i-th particle
 *  - Vector_t  specifying the electric field
 *  - Vector_t  specifying the magnetic field
 */

template <typename FieldFunction, typename ... Arguments>
class Stepper {
    
public:
    
    Stepper(const FieldFunction& fieldfunc) : fieldfunc_m(fieldfunc) { }
    
    virtual bool advance(PartBunchBase<double, 3>* bunch,
                         const size_t& i,
                         const double& t,
                         const double dt,
                         Arguments& ... args) const
    {
        bool isGood = doAdvance_m(bunch, i, t, dt, args...);

        bool isNaN = false;
        for (int j = 0; j < 3; ++j) {
            if (std::isnan(bunch->R[i](j)) ||
                std::isnan(bunch->P[i](j)) ||
                std::abs(bunch->R[i](j)) > 1.0e10 ||
                std::abs(bunch->P[i](j)) > 1.0e10) {
                isNaN = true;
                break;
            }
        }

        bool isBad = (!isGood || isNaN);
        if ( isBad ) {
            bunch->Bin[i] = -1;
        }
        return isBad;
    };
    virtual ~Stepper() {};

protected:
    const FieldFunction& fieldfunc_m;

private:
    virtual bool doAdvance_m(PartBunchBase<double, 3>* bunch,
                             const size_t& i,
                             const double& t,
                             const double dt,
                             Arguments& ... args) const = 0;
};

#endif
