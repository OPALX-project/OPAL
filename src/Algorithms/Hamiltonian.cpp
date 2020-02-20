#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(int truncOrder) : truncOrder_m(truncOrder)
{

    series_t::setGlobalTruncOrder(truncOrder+1);

    x     = series_t::makeVariable(0);
    px    = series_t::makeVariable(1);
    y     = series_t::makeVariable(2);
    py    = series_t::makeVariable(3);
    z     = series_t::makeVariable(4);
    delta = series_t::makeVariable(5);
}

Hamiltonian::series_t Hamiltonian::drift(const double& gamma0)
{
    double beta0 = this->getBeta_m(gamma0);

    return ( delta / beta0 )
            - sqrt((1./ beta0 + delta ) *(1./ beta0 + delta )
                    - ( px*px )
                    - ( py*py )
                    - 1./( beta0 * beta0 * gamma0 * gamma0 ),truncOrder_m+1);
}


Hamiltonian::series_t Hamiltonian::rbend(double& beta0,
                            double& gamma0,
                            double& /*q*/,
                            double& h,
                            double& k0)
{



        return ( delta / beta0 )
        - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                    - ( px*px )
                    - ( py*py )
                    - 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m+1
            ))
    - (h * x)
    * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                    - ( px*px )
                    - ( py*py )
                    - 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m
            ))
    + k0 * x * (1. + 0.5 * h* x);


}


Hamiltonian::series_t Hamiltonian::sbend(const double& gamma0,
                                         const double& h,
                                         const double& K0)
{
    double beta0 = this->getBeta_m(gamma0);

    return ( delta / beta0 )
                    - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                    - ( px*px )
                                    - ( py*py )
                                    - 1./( beta0*beta0 * gamma0*gamma0 ),(truncOrder_m+1)
                            ))
                    - (h * x)
                    * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                    - ( px*px )
                                    - ( py*py )
                                    - 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m
                            ))
                    + K0 * x * (1. + 0.5 * h* x);
}


Hamiltonian::series_t Hamiltonian::bendFringe(
                double& beta0,
                double& gamma0,
                double& /*h*/,
                double& /*k0*/,
                series_t& ax,
                series_t& az)
{
    if (truncOrder_m == 2){
        return ( delta / beta0 )
                        - sqrt((1./ beta0 + delta ) *(1./ beta0 + delta )
                                - ( px*px)
                                - ( py*py )
                                - 1./( beta0 * beta0 * gamma0 * gamma0 ),truncOrder_m+1
                        ) - az;
    }else{
        return ( delta / beta0 )
            - sqrt((1./ beta0 + delta ) *(1./ beta0 + delta )
                    - ( px*px - 2.0*px*ax - ax*ax)
                    - ( py*py )
                    - 1./( beta0 * beta0 * gamma0 * gamma0 ),truncOrder_m+1
            ) - az;
    }
    //std::cout << H << std::endl;
    //std::cout << H.getMaxOrder() << std::endl;
    /*H=( delta / beta0 )
                        - (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                        - ( px - ax )*( px - ax )
                                        - ( py*py )
                                        - 1./( beta0*beta0 * gamma0*gamma0 ),(truncOrder_m+1)))
                        - (h * x)
                        * (sqrt ((1./ beta0 + delta) *(1./ beta0 + delta)
                                        - ( px - ax )*( px - ax )
                                        - ( py*py )
                                        - 1./( beta0*beta0 * gamma0*gamma0 ), truncOrder_m))
                         + h * x * az;*/
}


Hamiltonian::series_t Hamiltonian::quadrupole(const double& gamma0,
                                              const double& /*q*/,
                                              const double& k1)
{
    double beta0 = this->getBeta_m(gamma0);

    return ( delta / beta0 )
    - sqrt ((1./ beta0 + delta ) *(1./ beta0 + delta)
            - ( px*px )
            - ( py*py )
            - 1./( beta0*beta0 * gamma0*gamma0 ),truncOrder_m+1
    )
    + 0.5 * k1 * (x*x - y*y);
}


Hamiltonian::series_t Hamiltonian::fringeField(const double& phi,
                                               const double& k0)
{
    if ( truncOrder_m > 1 ) {
        //TODO higher order thin lens approximation
    }
    // else

    double angle = phi;
    if ( k0 < 0.0 )
        angle = -phi;
    return -0.5 * (x * x - y * y) * k0 * std::tan(angle);
}


double Hamiltonian::getBeta_m(const double& gamma) {
    return std::sqrt(1.0 - 1.0 / (gamma * gamma) );
}
