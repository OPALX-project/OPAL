// -*- C++ -*-
/***************************************************************************
 *
 *
 * P3MPoissonSolver.hh
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

#ifndef P3M_POISSON_SOLVER_H_
#define P3M_POISSON_SOLVER_H_
const unsigned Dim = 3;

#ifdef dontOPTIMIZE_FIELD_ASSIGNMENT
#define FIELDASSIGNOPTIMIZATION __attribute__((optimize(0)))
#else
#define FIELDASSIGNOPTIMIZATION
#endif

#include <memory>
//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"
#include "Algorithms/PartBunchBase.h"

template <class T, unsigned Dim>
class PartBunchBase;

//////////////////////////////////////////////////////////////

class P3MPoissonSolver : public PoissonSolver {
public:
    // constructor and destructor
    P3MPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, double interaction_radius, double alpha, double eps);

    ~P3MPoissonSolver();

    void initFields();

    void calculateGridForces(PartBunchBase<double, 3> *bunch, double interaction_radius, double alpha, double eps);

    void calculatePairForces(PartBunchBase<double, 3> *bunch, double interaction_radius, double alpha, double eps);

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential with image charges at  -z
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    void computePotential(Field_t &rho, Vector_t hr);

    void applyConstantFocusing(PartBunchBase<double, 3> *bunch, double f, double r);
    void test(PartBunchBase<double, 3> *bunch);

    double getXRangeMin(unsigned short /*level*/) {return 1.0;}
    double getXRangeMax(unsigned short /*level*/) {return 1.0;}
    double getYRangeMin(unsigned short /*level*/) {return 1.0;}
    double getYRangeMax(unsigned short /*level*/) {return 1.0;}
    double getZRangeMin(unsigned short /*level*/) {return 1.0;}
    double getZRangeMax(unsigned short /*level*/) {return 1.0;}

    void computeAvgSpaceChargeForces(PartBunchBase<double, 3> *bunch);
    void compute_temperature(PartBunchBase<double, 3> *bunch);
    Inform &print(Inform &os) const;
private:

    BConds<double, Dim, Mesh_t, Center_t> bc_m;
    BConds<double, Dim, Mesh_t, Center_t> bcp_m;
    BConds<Vector_t, Dim, Mesh_t, Center_t> vbc_m;

    // rho_m is the charge-density field with mesh doubled in each dimension
    Field_t rho_m;
    Field_t phi_m;

    VField_t eg_m;

    // real field with layout of complex field: domain3_m
    Field_t greentr_m;

    CxField_t rhocmpl_m;
    CxField_t grncmpl_m;

    // grntr_m is the Fourier transformed Green's function
    // domain3_m and mesh3_ are used
    CxField_t grntr_m;

    // the FFT object
    std::unique_ptr<FFTC_t> fft_m;


    // Fields used to eliminate excess calculation in greensFunction()
    // mesh2_m and layout2_m are used
    IField_t grnIField_m[3];


    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;


    // tmp
    Field_t tmpgreen;

    // domains for the various fields
    NDIndex<3> domain_m;             // original domain, gridsize
    // mesh and gridsize defined outside of P3M class, given as


    NDIndex<3> domainP3MConstruct_m;


    double interaction_radius_m;
    double alpha_m;
    double eps_m;

    Vector_t hr_m;
    Vektor<int, 3> nr_m;

    // for tests
    Vektor<double,Dim> avgEF_m;
    double globSumEf_m[Dim];


public:
    Vektor<double,3> extend_l;
    Vektor<double,3> extend_r;


};

inline Inform &operator<<(Inform &os, const P3MPoissonSolver &fs) {
    return fs.print(os);
}



#endif

/***************************************************************************
 * $RCSfile: P3MPoissonSolver.hh,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2001/08/08 11:21:48 $
 ***************************************************************************/
