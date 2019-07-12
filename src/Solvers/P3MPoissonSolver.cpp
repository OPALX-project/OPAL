// -*- C++ -*-
/***************************************************************************
 *
 *
 * P3MPoissonSolver.cc
 *
 *
 *
 *
 *
 *
 *
 ***************************************************************************/

// include files
#include "Solvers/P3MPoissonSolver.h"
#include "Algorithms/PartBunch.h"
#include "Particle/BoxParticleCachingPolicy.h"
#include "Particle/PairBuilder/HashPairBuilderPeriodic.h"
#include "Particle/PairBuilder/HashPairBuilderPeriodicParallel.h"
//#include "Particle/PairBuilder/HashPairBuilderPeriodicParallel_globCHaining.h"
#include "Particle/PairBuilder/PairConditions.h"
#include "Structure/DataSink.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include <fstream>
//////////////////////////////////////////////////////////////////////////////
// a little helper class to specialize the action of the Green's function
// calculation.  This should be specialized for each dimension
// to the proper action for computing the Green's function.  The first
// template parameter is the full type of the Field to compute, and the second
// is the dimension of the data, which should be specialized.

//const double ke=1./(4.*M_PI*8.8e-14);
const double ke=2.532638e8;

template<unsigned int Dim>
struct SpecializedGreensFunction { };

template<>
struct SpecializedGreensFunction<3> {
    template<class T, class FT, class FT2>
    static void calculate(Vektor<T, 3> &hrsq, FT &grn, FT2 *grnI, double alpha,double eps) {
        double r;
        NDIndex<3> elem0=NDIndex<3>(Index(0,0), Index(0,0),Index(0,0));
        grn = grnI[0] * hrsq[0] + grnI[1] * hrsq[1] + grnI[2] * hrsq[2];
        NDIndex<3> lDomain_m = grn.getLayout().getLocalNDIndex();
        NDIndex<3> elem;
        for (int i=lDomain_m[0].min(); i<=lDomain_m[0].max(); ++i) {
            elem[0]=Index(i,i);
            for (int j=lDomain_m[1].min(); j<=lDomain_m[1].max(); ++j) {
                elem[1]=Index(j,j);
                for (int k=lDomain_m[2].min(); k<=lDomain_m[2].max(); ++k) {
                    elem[2]=Index(k,k);
                    r = real(sqrt(grn.localElement(elem)));
                    if(elem==elem0) {
                        //grn.localElement(elem) = ke*dcomplex(2*alpha/sqrt(M_PI)) ;
                        grn.localElement(elem) = 0 ;
                    }
                    else
                        grn.localElement(elem) = ke*dcomplex(erf(alpha*r)/(r+eps));
                }
            }
        }
    }
};

template<class T>
struct ApplyField {
    ApplyField(T c, double r, double epsilon, double alpha) : C(c), R(r), eps(epsilon), a(alpha) {}
    void operator()(std::size_t i, std::size_t j, PartBunch &P,Vektor<double,3> &shift) const
    {
        Vector_t diff = P.R[i] - (P.R[j]+shift);
        double sqr = 0;

        for (unsigned d = 0; d<Dim; ++d)
            sqr += diff[d]*diff[d];

        //compute r with softening parameter, unsoftened r is obtained by sqrt(sqr)
        if(sqr!=0) {
            double r = std::sqrt(sqr+eps*eps);
            //for order two transition
            if (P.Q[i]!=0 && P.Q[j]!=0) {
                //compute potential energy
                double phi =ke*(1.-erf(a*sqrt(sqr)))/r;

                //compute force
                Vector_t Fij = ke*C*(diff/sqrt(sqr))*((2.*a*exp(-a*a*sqr))/(sqrt(M_PI)*r)+(1.-erf(a*sqrt(sqr)))/(r*r));

                //Actual Force is F_ij multiplied by Qi*Qj
                //The electrical field on particle i is E=F/q_i and hence:
                P.Ef[i] -= P.Q[j]*Fij;
                P.Ef[j] += P.Q[i]*Fij;
                //update potential per particle
                P.Phi[i] += P.Q[j]*phi;
                P.Phi[j] += P.Q[i]*phi;
            }
        }
    }
    T C;
    double R;
    double eps;
    double a;
};


////////////////////////////////////////////////////////////////////////////

// constructor


P3MPoissonSolver::P3MPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, double interaction_radius, double alpha, double eps):
    mesh_m(mesh),
    layout_m(fl),
    interaction_radius_m(interaction_radius),
    alpha_m(alpha),
    eps_m(eps)
{
    Inform msg("P3MPoissonSolver::P3MPoissonSolver ");

    domain_m = layout_m->getDomain();

    for (unsigned int i = 0; i < 2 * Dim; ++i) {
        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new ParallelPeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new PeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
    }

    for (unsigned int i = 0; i < 2 * Dim; ++i) {
        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new ParallelPeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new PeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
    }


    // For efficiency in the FFT's, we can use a parallel decomposition
    // which can be serial in the first dimension.
    // e_dim_tag decomp[3];
    // for(int d = 0; d < 3; ++d) {
    //     decomp[d] = layout_m->getRequestedDistribution(d);
    // }

    // The FFT's require double-sized field sizes in order to (more closely
    // do not understand this ...)
    // simulate an isolated system.  The FFT of the charge density field, rho,
    // would otherwise mimic periodic boundary conditions, i.e. as if there were
    // several beams set a periodic distance apart.  The double-sized fields
    // alleviate this problem.
    for(unsigned int i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
    }

    rhocmpl_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));
    grncmpl_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));

    rho_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1), bc_m);
    phi_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1), bcp_m);
    eg_m.initialize (*mesh_m, *layout_m, GuardCellSizes<Dim>(1), vbc_m);

    bool compressTemps = true;
    // create the FFT object
    fft_m = std::unique_ptr<FFTC_t>(new FFTC_t(layout_m->getDomain(), compressTemps));
    fft_m->setDirectionName(+1, "forward");
    fft_m->setDirectionName(-1, "inverse");
}

void P3MPoissonSolver::initFields() {

    domain_m = layout_m->getDomain();

    for (unsigned int i = 0; i < 2 * Dim; ++i) {
        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new ParallelPeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new PeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
    }

    for (unsigned int i = 0; i < 2 * Dim; ++i) {
        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new ParallelPeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new PeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
    }

    for(unsigned int i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
    }

    rhocmpl_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));
    grncmpl_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));

    rho_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1), bc_m);
    phi_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1), bcp_m);
    eg_m.initialize (*mesh_m, *layout_m, GuardCellSizes<Dim>(1), vbc_m);

    bool compressTemps = true;
    if (fft_m)
        fft_m.reset();

    // create the FFT object
    fft_m = std::unique_ptr<FFTC_t>(new FFTC_t(layout_m->getDomain(), compressTemps));
    fft_m->setDirectionName(+1, "forward");
    fft_m->setDirectionName(-1, "inverse");

}



////////////////////////////////////////////////////////////////////////////
// destructor
P3MPoissonSolver::~P3MPoissonSolver() {
}

void P3MPoissonSolver::calculatePairForces(PartBunchBase<double, 3> *bunch, double interaction_radius, double alpha, double eps) {
    if (interaction_radius>0){
        if (Ippl::getNodes() > 1) {
            PartBunch &tmpBunch = *(dynamic_cast<PartBunch*>(bunch));
            HashPairBuilderPeriodicParallel<PartBunch> HPB(tmpBunch);
            HPB.for_each(RadiusCondition<double, Dim>(interaction_radius), ApplyField<double>(-1,interaction_radius,eps,alpha),extend_l, extend_r);
        }
        else {
            PartBunch &tmpBunch = *(dynamic_cast<PartBunch*>(bunch));
            HashPairBuilderPeriodic<PartBunch> HPB(tmpBunch);
            HPB.for_each(RadiusCondition<double, Dim>(interaction_radius), ApplyField<double>(-1,interaction_radius,eps,alpha),extend_l, extend_r);
        }
    }

}

void P3MPoissonSolver::calculateGridForces(PartBunchBase<double, 3> *bunch, double interaction_radius, double alpha, double eps){

    Inform msg ("calculateGridForces ");
    Vector_t l,h;

    // (1) scatter charge to charge density grid and transform to fourier space
    rho_m[domain_m]=0;
    bunch->Q.scatter(rho_m, bunch->R, IntrplCIC_t());

    rhocmpl_m[domain_m] = rho_m[domain_m]/(hr_m[0]*hr_m[1]*hr_m[2]);

    fft_m->transform("inverse", rhocmpl_m);

    // (2) compute Greens function in real space and transform to fourier space
    IField_t grnIField[3];

    // This loop stores in grnIField_m[i] the index of the ith dimension mirrored at the central axis. e.g.
    // grnIField_m[0]=[(0 1 2 3 ... 3 2 1) ; (0 1 2 3 ... 3 2 1; ...)]
    for (int i = 0; i < 3; ++i) {
        grnIField[i].initialize(*mesh_m, *layout_m);
        grnIField[i][domain_m] = where(lt(domain_m[i], nr_m[i]/2),
                                       domain_m[i] * domain_m[i],
                                       (nr_m[i]-domain_m[i]) *
                                       (nr_m[i]-domain_m[i]));
    }
    Vector_t hrsq(hr_m * hr_m);
    SpecializedGreensFunction<3>::calculate(hrsq, grncmpl_m, grnIField, alpha, eps);

    //transform G -> Ghat and store in grncmpl_m
    fft_m->transform("inverse", grncmpl_m);
    //multiply in fourier space and obtain PhiHat in rhocmpl_m
    rhocmpl_m *= grncmpl_m;

    // (3) Backtransformation: compute potential field in real space and E=-Grad Phi
    //compute electrostatic potential Phi in real space by FFT PhiHat -> Phi and store it in rhocmpl_m
    fft_m->transform("forward", rhocmpl_m);

    //take only the real part and store in phi_m (has periodic bc instead of interpolation bc)
    phi_m = real(rhocmpl_m)*hr_m[0]*hr_m[1]*hr_m[2];
    //dumpVTKScalar( phi_m, this,it, "Phi_m") ;

    //compute Electric field on the grid by -Grad(Phi) store in eg_m
    eg_m = -Grad1Ord(phi_m, eg_m);

    //interpolate the electric field to the particle positions
    bunch->Ef.gather(eg_m, bunch->R,  IntrplCIC_t());
    //interpolate electrostatic potenital to the particle positions
    bunch->Phi.gather(phi_m, bunch->R, IntrplCIC_t());
}



////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric potential from the image charge by solving
// the Poisson's equation

void P3MPoissonSolver::computePotential(Field_t &rho, Vector_t hr, double zshift) {


}

void P3MPoissonSolver::computeAvgSpaceChargeForces(PartBunchBase<double, 3> *bunch) {

    Inform m("computeAvgSpaceChargeForces ");

    const double N =  static_cast<double>(bunch->getTotalNum());
    double locAvgEf[Dim]={};
    for (unsigned i=0; i<bunch->getLocalNum(); ++i) {
        locAvgEf[0]+=fabs(bunch->Ef[i](0));
        locAvgEf[1]+=fabs(bunch->Ef[i](1));
        locAvgEf[2]+=fabs(bunch->Ef[i](2));
    }

    reduce(&(locAvgEf[0]), &(locAvgEf[0]) + Dim,
           &(globSumEf_m[0]), OpAddAssign());

    //    m << "globSumEF = " << globSumEf_m[0] << "\t" << globSumEf_m[1] << "\t" << globSumEf_m[2] << endl;

    avgEF_m[0]=globSumEf_m[0]/N;
    avgEF_m[1]=globSumEf_m[1]/N;
    avgEF_m[2]=globSumEf_m[2]/N;

}


void P3MPoissonSolver::applyConstantFocusing(PartBunchBase<double, 3> *bunch, double f, double r){
    Vektor<double,Dim> beam_center(0,0,0);
    Vector_t Rrel;
    double scFoc = sqrt(dot(avgEF_m,avgEF_m));
    for (unsigned i=0; i<bunch->getLocalNum(); ++i) {
        Rrel=bunch->R[i] - beam_center;
        bunch->Ef[i] += Rrel/r*f*scFoc;
    }

}

// given a charge-density field rho and a set of mesh spacings hr,
// compute the scalar potential in open space
void P3MPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {


}

void P3MPoissonSolver::compute_temperature(PartBunchBase<double, 3> *bunch) {

    Inform m("compute_temperature ");
    // double loc_temp[Dim]={0.0,0.0,0.0};
    double loc_avg_vel[Dim]={0.0,0.0,0.0};
    double avg_vel[Dim]={0.0,0.0,0.0};

    for(unsigned long k = 0; k < bunch->getLocalNum(); ++k) {
        for(unsigned i = 0; i < Dim; i++) {
            loc_avg_vel[i]   += bunch->P[k](i);
        }
    }
    reduce(&(loc_avg_vel[0]), &(loc_avg_vel[0]) + Dim,
           &(avg_vel[0]), OpAddAssign());

    const double N =  static_cast<double>(bunch->getTotalNum());
    avg_vel[0]=avg_vel[0]/N;
    avg_vel[1]=avg_vel[1]/N;
    avg_vel[2]=avg_vel[2]/N;

    m << "avg_vel[0]= " << avg_vel[0] << " avg_vel[1]= " << avg_vel[1]  << " avg_vel[2]= " << avg_vel[2] << endl;

    /*
      for(unsigned long k = 0; k < bunch.getLocalNum(); ++k) {
      for(unsigned i = 0; i < Dim; i++) {
      loc_temp[i]   += (bunch.P[k](i)-avg_vel[i])*(bunch.P[k](i)-avg_vel[i]);
      }
      }
      reduce(&(loc_temp[0]), &(loc_temp[0]) + Dim,
      &(temperature[0]), OpAddAssign());
      temperature[0]=temperature[0]/N;
      temperature[1]=temperature[1]/N;
      temperature[2]=temperature[2]/N;
    */
}

void P3MPoissonSolver::test(PartBunchBase<double, 3> *bunch) {
    Inform msg("P3MPoissonSolver::test ");

    // set special conditions for this test
    const double mi = 1.0;
    const double qi = -1.0;
    const double qom = qi/mi;
    const double beam_radius = 0.001774;
    const double f = 1.5;
    const double dt = bunch->getdT();

    OpalData *opal = OpalData::getInstance();
    DataSink *ds = opal->getDataSink();

    //    std::vector<std::pair<std::string, unsigned int> > collimatorLosses; // just empty
    Vector_t FDext[6];

    bunch->Q = qi;
    bunch->M = mi;

    bunch->calcBeamParameters();

    initFields();

    for (int i=0; i<3; i++) {
        extend_r[i] =  hr_m[i]*nr_m[i]/2;
        extend_l[i] = -hr_m[i]*nr_m[i]/2;
    }

    msg << *this << endl;

    // calculate initial space charge forces
    calculateGridForces(bunch, interaction_radius_m, alpha_m, eps_m);
    calculatePairForces(bunch, interaction_radius_m, alpha_m, eps_m);

    //avg space charge forces for constant focusing
    computeAvgSpaceChargeForces(bunch);

    for (int it=0; it<1000; it++) {

        // advance the particle positions
        // basic leapfrogging timestep scheme.  velocities are offset
        // by half a timestep from the positions.

        assign(bunch->R, bunch->R + dt * bunch->P);

        bunch->update();

        calculateGridForces(bunch, interaction_radius_m, alpha_m, eps_m);
        calculatePairForces(bunch, interaction_radius_m, alpha_m, eps_m);
        applyConstantFocusing(bunch,f,beam_radius);

        assign(bunch->P, bunch->P + dt * qom * bunch->Ef);

        if (it%10 == 0){
            bunch->calcBeamParameters();
            ds->dumpSDDS(bunch, FDext, it);
        }
        msg << "Finished iteration " << it << endl;
    }
}


Inform &P3MPoissonSolver::print(Inform &os) const {
    os << "* ************* P 3 M - P o i s s o n S o l v e r *************** " << endl;
    os << "* h        " << hr_m << '\n';
    os << "* RC       " << interaction_radius_m << '\n';
    os << "* ALPHA    " << alpha_m << '\n';
    os << "* EPSILON  " << eps_m << '\n';
    os << "* Extend L " << extend_l << '\n';
    os << "* Extend R " << extend_r << '\n';
    os << "* hr       " << hr_m << '\n';
    os << "* nr       " << nr_m << '\n';
    os << "* *************************************************************** " << endl;
    return os;
}

/***************************************************************************
 * $RCSfile: P3MPoissonSolver.cc,v $   $Author: adelmann $
 * $Revision: 1.6 $   $Date: 2001/08/16 09:36:08 $
 ***************************************************************************/
