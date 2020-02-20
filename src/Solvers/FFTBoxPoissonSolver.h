// -*- C++ -*-
/***************************************************************************
 *
 *
 * FFTBoxPoissonSolver.hh
 *
 * Open BC in x,y and z.
 *
 *
 *
 *
 *
 *
 ***************************************************************************/

////////////////////////////////////////////////////////////////////////////
// This class contains methods for solving Poisson's equation for the
// space charge portion of the calculation.
////////////////////////////////////////////////////////////////////////////

#ifndef FFT_BOXPOISSON_SOLVER_H_
#define FFT_BOXPOISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"
class PartBunch;
//////////////////////////////////////////////////////////////

class FFTBoxPoissonSolver : public PoissonSolver {
public:
    // constructor and destructor
    FFTBoxPoissonSolver(PartBunch &bunch, std::string greensFuntion);

    FFTBoxPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::string greensFunction, double boxSize);

    ~FFTBoxPoissonSolver();

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential with image charges at  -z
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    void computePotential(Field_t &rho, Vector_t hr);

    // compute the green's function for a Poisson problem and put it in in grntm_m
    // uses grnIField_m to eliminate excess calculation in greenFunction()
    // given mesh information in nr and hr
    void greensFunction();

    /// compute the integrated Green function as described in <A HREF="http://prst-ab.aps.org/abstract/PRSTAB/v9/i4/e044204">Three-dimensional quasistatic model for high brightness beam dynamics simulation</A> by Qiang et al.
    void integratedGreensFunction();

    /// compute the shifted integrated Green function as described in <A HREF="http://prst-ab.aps.org/abstract/PRSTAB/v9/i4/e044204">Three-dimensional quasistatic model for high brightness beam dynamics simulation</A> by Qiang et al.
    void shiftedIntGreensFunction(double zshift);

    double getXRangeMin(unsigned short /*level*/) {return -a_m;}
    double getXRangeMax(unsigned short /*level*/) {return  a_m;}
    double getYRangeMin(unsigned short /*level*/) {return -a_m;}
    double getYRangeMax(unsigned short /*level*/) {return  a_m;}
    double getZRangeMin(unsigned short /*level*/) {return -a_m; }
    double getZRangeMax(unsigned short /*level*/) {return  a_m; }
        void test(PartBunchBase<double, 3> */*bunch*/) { }


    Inform &print(Inform &os) const;

private:

    // rho2_m is the charge-density field with mesh doubled in each dimension
    Field_t rho2_m;

    // grntr_m is the Fourier transformed Green's function
    Field_t grntr_m;

    // Fields used to eliminate excess calculation in greensFunction()
    // mesh2_m and layout2_m are used
    IField_t grnIField_m[3];

    SINE_t *sine_m;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;

    // mesh and layout objects for rho2_m
    Mesh_t *mesh2_m;
    FieldLayout_t *layout2_m;

    // tmp
    Field_t tmpgreen;

    // domains for the various fields
    NDIndex<3> domain_m;             // original domain, gridsize
    // mesh and gridsize defined outside of FFT class, given as
    // parameter to the constructor (mesh and layout object).
    NDIndex<3> domain2_m;            // doubled gridsize (2*Nx,2*Ny,2*Nz)
    // (2*Nx,Ny,2*Nz)
    // mesh spacing and size values
    Vector_t hr_m;
    Vektor<int, 3> nr_m;

    std::string greensFunction_m;

    double a_m; // the box size

    IpplTimings::TimerRef GreensFunctionTimer_m;

    IpplTimings::TimerRef IntGreensFunctionTimer1_m;
    IpplTimings::TimerRef IntGreensFunctionTimer2_m;
    IpplTimings::TimerRef IntGreensFunctionTimer3_m;
    IpplTimings::TimerRef IntGreensFunctionTimer4_m;

    IpplTimings::TimerRef ShIntGreensFunctionTimer1_m;
    IpplTimings::TimerRef ShIntGreensFunctionTimer2_m;
    IpplTimings::TimerRef ShIntGreensFunctionTimer3_m;
    IpplTimings::TimerRef ShIntGreensFunctionTimer4_m;

    IpplTimings::TimerRef GreensFunctionTimer1_m;
    IpplTimings::TimerRef GreensFunctionTimer4_m;
};

inline Inform &operator<<(Inform &os, const FFTBoxPoissonSolver &fs) {
    return fs.print(os);
}



#endif