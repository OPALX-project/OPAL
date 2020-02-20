// -*- C++ -*-
/***************************************************************************
 *
 *
 * FFTPoissonSolver.hh
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

#ifndef FFT_POISSON_SOLVER_H_
#define FFT_POISSON_SOLVER_H_


#ifdef dontOPTIMIZE_FIELD_ASSIGNMENT
#define FIELDASSIGNOPTIMIZATION __attribute__((optimize(0)))
#else
#define FIELDASSIGNOPTIMIZATION
#endif

#include <memory>
//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"

#ifdef OPAL_DKS
#include "DKSOPAL.h"
#endif

class PartBunch;

//////////////////////////////////////////////////////////////

class FFTPoissonSolver : public PoissonSolver {
public:
    // constructor and destructor
    FFTPoissonSolver(PartBunch &bunch, std::string greensFuntion);

    FFTPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::string greensFunction, std::string bcz);

    ~FFTPoissonSolver();

#ifdef OPAL_DKS
    DKSOPAL dksbase;
#endif

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential with image charges at  -z
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    void computePotential(Field_t &rho, Vector_t hr);

    void computePotentialDKS(Field_t &rho);
    // compute the green's function for a Poisson problem and put it in in grntm_m
    // uses grnIField_m to eliminate excess calculation in greenFunction()
    // given mesh information in nr and hr
    void greensFunction();

    /// compute the integrated Green function as described in <A HREF="http://prst-ab.aps.org/abstract/PRSTAB/v9/i4/e044204">Three-dimensional quasistatic model for high brightness beam dynamics simulation</A> by Qiang et al.
    void integratedGreensFunction();

    /// Uses DKS to offload the computation of Greens function on the GPU
    //compute the integrated Green function as described in <A HREF="http://prst-ab.aps.org/abstract/PRSTAB/v9/i4/e044204">Three-dimensional quasistatic model for high brightness beam dynamics simulation</A> by Qiang et al.
    void integratedGreensFunctionDKS();

    /// compute the shifted integrated Green function as described in <A HREF="http://prst-ab.aps.org/abstract/PRSTAB/v9/i4/e044204">Three-dimensional quasistatic model for high brightness beam dynamics simulation</A> by Qiang et al.
    void shiftedIntGreensFunction(double zshift);

    double getXRangeMin(unsigned short /*level*/) {return 1.0;}
    double getXRangeMax(unsigned short /*level*/) {return 1.0;}
    double getYRangeMin(unsigned short /*level*/) {return 1.0;}
    double getYRangeMax(unsigned short /*level*/) {return 1.0;}
    double getZRangeMin(unsigned short /*level*/) {return 1.0;}
    double getZRangeMax(unsigned short /*level*/) {return 1.0;}
    void test(PartBunchBase<double, 3> */*bunch*/) { }

    Inform &print(Inform &os) const;
private:
    void initializeFields();
    void mirrorRhoField() FIELDASSIGNOPTIMIZATION;
    void mirrorRhoField(Field_t & ggrn2);// FIELDASSIGNOPTIMIZATION;

    // rho2_m is the charge-density field with mesh doubled in each dimension
    Field_t rho2_m;

    // real field with layout of complex field: domain3_m
    Field_t greentr_m;

    // rho2tr_m is the Fourier transformed charge-density field
    // domain3_m and mesh3_ are used
    CxField_t rho2tr_m;
    CxField_t imgrho2tr_m;

    // grntr_m is the Fourier transformed Green's function
    // domain3_m and mesh3_ are used
    CxField_t grntr_m;

#ifdef OPAL_DKS
    //pointer for Fourier transformed Green's function on GPU
    void * grntr_m_ptr;
    void * rho2_m_ptr;
    void * tmpgreen_ptr;
    void *rho2real_m_ptr;
    void *rho2tr_m_ptr;

    //stream id for calculating greens function
    int streamGreens;
    int streamFFT;

#endif

    // Fields used to eliminate excess calculation in greensFunction()
    // mesh2_m and layout2_m are used
    IField_t grnIField_m[3];

    // the FFT object
    std::unique_ptr<FFT_t> fft_m;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;

    // mesh and layout objects for rho2_m
    std::unique_ptr<Mesh_t> mesh2_m;
    std::unique_ptr<FieldLayout_t> layout2_m;

    //
    std::unique_ptr<Mesh_t> mesh3_m;
    std::unique_ptr<FieldLayout_t> layout3_m;

    // mesh and layout for integrated greens function
    std::unique_ptr<Mesh_t> mesh4_m;
    std::unique_ptr<FieldLayout_t> layout4_m;

    // tmp
    Field_t tmpgreen_m;

    // domains for the various fields
    NDIndex<3> domain_m;             // original domain, gridsize
    // mesh and gridsize defined outside of FFT class, given as
    // parameter to the constructor (mesh and layout object).
    NDIndex<3> domain2_m;            // doubled gridsize (2*Nx,2*Ny,2*Nz)
    NDIndex<3> domain3_m;            // field for the complex values of the RC transformation
    NDIndex<3> domain4_m;
    // (2*Nx,Ny,2*Nz)
    NDIndex<3> domainFFTConstruct_m;

    // mesh spacing and size values
    Vector_t hr_m;
    Vektor<int, 3> nr_m;

    /// for defining the boundary conditions
    BConds<double, 3, Mesh_t, Center_t> bc_m;
    BConds<Vector_t, 3, Mesh_t, Center_t> vbc_m;

    bool bcz_m;
    bool integratedGreens_m;
    IpplTimings::TimerRef GreensFunctionTimer_m;
  /*
    IpplTimings::TimerRef IntGreensFunctionTimer1_m;
    IpplTimings::TimerRef IntGreensFunctionTimer2_m;
    IpplTimings::TimerRef IntGreensFunctionTimer3_m;
    IpplTimings::TimerRef IntGreensMirrorTimer1_m;

    IpplTimings::TimerRef ShIntGreensFunctionTimer1_m;
    IpplTimings::TimerRef ShIntGreensFunctionTimer2_m;
    IpplTimings::TimerRef ShIntGreensFunctionTimer3_m;
    IpplTimings::TimerRef ShIntGreensFunctionTimer4_m;
    IpplTimings::TimerRef IntGreensMirrorTimer2_m;

    IpplTimings::TimerRef GreensFunctionTimer1_m;
    IpplTimings::TimerRef GreensFunctionTimer4_m;
  */
    IpplTimings::TimerRef ComputePotential_m;
};

inline Inform &operator<<(Inform &os, const FFTPoissonSolver &fs) {
    return fs.print(os);
}



#endif