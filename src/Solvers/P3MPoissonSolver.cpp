//
// Class P3MPoissonSolver
//   This class contains methods for solving Poisson's equation for the
//   space charge portion of the calculation including collisions.
//
// Copyright (c) 2016, Benjamin Ulmer, ETH Zürich
//               2022, Sriramkrishnan Muralikrishnan, PSI
// All rights reserved
//
// Implemented as part of the Master thesis
// "The P3M Model on Emerging Computer Architectures With Application to Microbunching"
// (http://amas.web.psi.ch/people/aadelmann/ETH-Accel-Lecture-1/projectscompleted/cse/thesisBUlmer.pdf)
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

#include "Solvers/P3MPoissonSolver.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunch.h"
#include "Particle/BoxParticleCachingPolicy.h"
#include "Particle/PairBuilder/HashPairBuilderParallel.h"
#include "Particle/PairBuilder/PairConditions.h"
#include "Physics/Physics.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <fstream>

//////////////////////////////////////////////////////////////////////////////
// a little helper class to specialize the action of the Green's function
// calculation.  This should be specialized for each dimension
// to the proper action for computing the Green's function.  The first
// template parameter is the full type of the Field to compute, and the second
// is the dimension of the data, which should be specialized.


template<unsigned int Dim>
struct P3MGreensFunction { };

template<>
struct P3MGreensFunction<3> {
    template<class T, class FT, class FT2>
    static void calculate(Vektor<T, 3>& hrsq_r, FT& grn_r, FT2* grnI_p, double alpha) {
        grn_r = grnI_p[0] * hrsq_r[0] + grnI_p[1] * hrsq_r[1] + grnI_p[2] * hrsq_r[2];
        grn_r = erf(alpha*sqrt(grn_r))/(sqrt(grn_r));
        grn_r[0][0][0] = grn_r[0][0][1];
    }
};

template<class T>
struct ApplyField {
    ApplyField(double alpha_, double ke_, bool isIntGreen_) : alpha(alpha_), ke(ke_), 
                                                             isIntGreen(isIntGreen_) {}
    void operator()(std::size_t i, std::size_t j, PartBunch& P_r) const
    {
        Vector_t diff = P_r.R[i] - P_r.R[j];
        double sqr = 0;

        for (unsigned d = 0; d<Dim; ++d) {
            sqr += diff[d]*diff[d];
        }

        if(sqr!=0) {
            double r = std::sqrt(sqr);

            //compute force
            Vector_t Fij;

            //Differentiate the PP Green's function 
            //(1-erf(\alpha r))/r (for standard)
            //(1/(2r)) * (\xi^3 - 3\xi +2) (for integrated) where 
            //xi = r/interaction_radius and multiply it with r unit vector
            double xi = r/alpha;
            Fij = ((double)isIntGreen * (-ke*(diff/(2*sqr))*((std::pow(xi,3) - 3*xi + 2)/r 
                   + 3*(1 - std::pow(xi,2))/alpha)))
                   + ((double)(1.0 - isIntGreen) 
                   * (-ke*(diff/r)*((2.*alpha*std::exp(-alpha*alpha*sqr))
                   / (std::sqrt(Physics::pi)*r) + (1.-std::erf(alpha*r))/(r*r))));

            //Actual Force is F_ij multiplied by Qi*Qj
            //The electrical field on particle i is E=F/q_i and hence:
            P_r.Ef[i] -= P_r.Q[j]*Fij;
            P_r.Ef[j] += P_r.Q[i]*Fij;
        }
    
    }
    double alpha;
    double ke;
    bool isIntGreen;
};


////////////////////////////////////////////////////////////////////////////

// constructor

P3MPoissonSolver::P3MPoissonSolver(Mesh_t* mesh_p, FieldLayout_t* fl_p, 
                                   double interaction_radius, 
                                   double alpha, std::string greensFunction):
    mesh_mp(mesh_p),
    layout_mp(fl_p),
    interaction_radius_m(interaction_radius),
    alpha_m(alpha)
{
    Inform msg("P3MPoissonSolver::P3MPoissonSolver ");

    integratedGreens_m = (greensFunction == std::string("INTEGRATED"));
    initializeFields();
    ke_m = 1.0 / (4 * Physics::pi * Physics::epsilon_0);

    GreensFunctionTimer_m = IpplTimings::getTimer("GreensFTotalP3M");
    ComputePotential_m = IpplTimings::getTimer("ComputePotentialP3M");
    CalculatePairForces_m = IpplTimings::getTimer("CalculatePairForcesP3M");
}

////////////////////////////////////////////////////////////////////////////
// destructor
P3MPoissonSolver::~P3MPoissonSolver() {
}


void P3MPoissonSolver::initializeFields() {

    domain_m = layout_mp->getDomain();

    // For efficiency in the FFT's, we can use a parallel decomposition
    // which can be serial in the first dimension.
    e_dim_tag decomp[3];
    e_dim_tag decomp2[3];
    for(int d = 0; d < 3; ++ d) {
        decomp[d] = layout_mp->getRequestedDistribution(d);
        decomp2[d] = layout_mp->getRequestedDistribution(d);
    }

    // The FFT's require double-sized field sizes in order to
    // simulate an isolated system.  The FFT of the charge density field, rho,
    // would otherwise mimic periodic boundary conditions, i.e. as if there were
    // several beams set a periodic distance apart.  The double-sized fields
    // alleviate this problem.
    for (int i = 0; i < 3; ++ i) {
        hr_m[i] = mesh_mp->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
        domain2_m[i] = Index(2 * nr_m[i] + 1);
    }

    // create double sized mesh and layout objects for the use in the FFT's
    mesh2_mp = std::unique_ptr<Mesh_t>(new Mesh_t(domain2_m));
    layout2_mp = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh2_mp, decomp));

    rho2_m.initialize(*mesh2_mp, *layout2_mp);

    // Create the domain for the transformed (complex) fields.  Do this by
    // taking the domain from the doubled mesh, permuting it to the right, and
    // setting the 2nd dimension to have n/2 + 1 elements.
    domain3_m[0] = Index(2 * nr_m[2] + 1);
    domain3_m[1] = Index(nr_m[0] + 2);
    domain3_m[2] = Index(2 * nr_m[1] + 1);

    // create mesh and layout for the new real-to-complex FFT's, for the
    // complex transformed fields
    mesh3_mp = std::unique_ptr<Mesh_t>(new Mesh_t(domain3_m));
    layout3_mp = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh3_mp, decomp2));

    rho2tr_m.initialize(*mesh3_mp, *layout3_mp);
    grntr_m.initialize(*mesh3_mp, *layout3_mp);

    for (int i = 0; i < 3; ++ i) {
        domain4_m[i] = Index(nr_m[i] + 2);
    }
    mesh4_mp = std::unique_ptr<Mesh_t>(new Mesh_t(domain4_m));
    layout4_mp = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh4_mp, decomp));

    tmpgreen_m.initialize(*mesh4_mp, *layout4_mp);

    // create a domain used to indicate to the FFT's how to construct it's
    // temporary fields.  This is the same as the complex field's domain,
    // but permuted back to the left.
    NDIndex<3> tmpdomain;
    tmpdomain = layout3_mp->getDomain();
    for (int i = 0; i < 3; ++ i)
        domainFFTConstruct_m[i] = tmpdomain[(i+1) % 3];

    // create the FFT object
    fft_mp = std::unique_ptr<FFT_t>(new FFT_t(layout2_mp->getDomain(), 
                                             domainFFTConstruct_m));

    if(!integratedGreens_m) {
        // these are fields that are used for calculating the Green's function.
        // they eliminate some calculation at each time-step.
        for (int i = 0; i < 3; ++ i) {
            grnIField_m[i].initialize(*mesh2_mp, *layout2_mp);
            grnIField_m[i][domain2_m] = where(lt(domain2_m[i], nr_m[i]),
                                              domain2_m[i] * domain2_m[i],
                                              (2 * nr_m[i] - domain2_m[i]) *
                                              (2 * nr_m[i] - domain2_m[i]));
        }
    }
}


void P3MPoissonSolver::calculatePairForces(PartBunchBase<double, 3>* bunch_p, double gammaz) {
    
    IpplTimings::startTimer(CalculatePairForces_m);
    if (interaction_radius_m>0){
        PartBunch& tmpBunch_r = *(dynamic_cast<PartBunch*>(bunch_p));
        std::size_t size = tmpBunch_r.getLocalNum()+tmpBunch_r.getGhostNum();
       
        //Take the particles to the boosted frame
        for(std::size_t i = 0;i<size;++i)
        {
            tmpBunch_r.R[i](2) = tmpBunch_r.R[i](2) * gammaz;
        }
        
        HashPairBuilderParallel<PartBunch> HPB(tmpBunch_r,gammaz);
        if(integratedGreens_m) {
            //Note: alpha_m is not used for the integrated Green's function
            //approach
            HPB.forEach(RadiusCondition<double, Dim>(interaction_radius_m), 
                         ApplyField<double>(interaction_radius_m,ke_m,integratedGreens_m));
        }
        else {
            HPB.forEach(RadiusCondition<double, Dim>(interaction_radius_m), 
                         ApplyField<double>(alpha_m,ke_m,integratedGreens_m));
        }
        
        //Bring the particles to the lab frame
        for(std::size_t i = 0;i<size;++i)
        {
            tmpBunch_r.R[i](2) = tmpBunch_r.R[i](2) / gammaz;
        }
    }
    IpplTimings::stopTimer(CalculatePairForces_m);
}


// given a charge-density field rho and a set of mesh spacings hr,
// compute the scalar potential in open space
void P3MPoissonSolver::computePotential(Field_t& rho_r, Vector_t hr) {
    IpplTimings::startTimer(ComputePotential_m);

    // use grid of complex doubled in both dimensions
    // and store rho in lower left quadrant of doubled grid
    rho2_m = 0.0;

    rho2_m[domain_m] = rho_r[domain_m];

    // needed in greens function
    hr_m = hr;

    // FFT double-sized charge density
    // we do a backward transformation so that we dont have to 
    // account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    fft_mp->transform(-1, rho2_m, rho2tr_m);

    // must be called if the mesh size has changed
    IpplTimings::startTimer(GreensFunctionTimer_m);
    if(integratedGreens_m)
        integratedGreensFunction();
    else
        greensFunction();
    IpplTimings::stopTimer(GreensFunctionTimer_m);
    // multiply transformed charge density
    // and transformed Green function
    // Don't divide by (2*nx_m)*(2*ny_m), as Ryne does;
    // this normalization is done in IPPL's fft routine.
    rho2tr_m *= grntr_m;

    // inverse FFT, rho2_m equals to the electrostatic potential
    fft_mp->transform(+1, rho2tr_m, rho2_m);
    // end convolution

    // back to physical grid
    // reuse the charge density field to store the electrostatic potential
    rho_r[domain_m] = rho2_m[domain_m];


    rho_r *= hr[0] * hr[1] * hr[2];

    IpplTimings::stopTimer(ComputePotential_m);
}

////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric potential from the image charge by solving
// the Poisson's equation

void P3MPoissonSolver::computePotential(Field_t& /*rho*/, Vector_t /*hr*/, double /*zshift*/) {

    throw OpalException("P3MPoissonSolver", "not implemented yet");

}

void P3MPoissonSolver::greensFunction() {

    Vector_t hrsq(hr_m * hr_m);
    P3MGreensFunction<3>::calculate(hrsq, rho2_m, grnIField_m, alpha_m);

    // we do a backward transformation so that we dont have to account 
    // for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    fft_mp->transform(-1, rho2_m, grntr_m);
}

/** If the beam has a longitudinal size >> transverse size the
 * direct Green function at each mesh point is not efficient
 * (needs a lot of mesh points along the transverse size to
 * get a good resolution)
 *
 * If we assume the charge density function is uniform within 
 * each cell then we can integrate the Green's function within
 * a cell and use it to improve the accuracy.
 *
 * We do not use the standard P3M Green's function
 * erf(\alpha r)/r and do the integration as no closed form
 * expression is available. Instead we use the zeroth order
 * truncated polynomial from "Hünenberger, P. H. (2000). 
 * The Journal of Chemical Physics, 113(23), 10464-10476" 
 * Table I in the appendix. The integration of higher-order 
 * polynomials gives complex Green's functions and hence 
 * has not been implemented yet.
 */
void P3MPoissonSolver::integratedGreensFunction() {
    /**
     * This integral can be calculated analytically in a closed from:
     */
    NDIndex<3> idx =  layout4_mp->getLocalNDIndex();
    double cellVolume = hr_m[0] * hr_m[1] * hr_m[2];
    tmpgreen_m = 0.0;

    for(int k = idx[2].first(); k <= idx[2].last() + 1; k++) {
        for(int j = idx[1].first(); j <=  idx[1].last() + 1; j++) {
            for(int i = idx[0].first(); i <= idx[0].last() + 1; i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = k * hr_m[2] - hr_m[2] / 2;

                double r = std::sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                double tmpgrn = 0.0;
                if(r >= interaction_radius_m) {
                    tmpgrn  = -vv(2) * vv(2) * std::atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                    tmpgrn += -vv(1) * vv(1) * std::atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                    tmpgrn += -vv(0) * vv(0) * std::atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                    tmpgrn += vv(1) * vv(2) * std::log(vv(0) + r);
                    tmpgrn += vv(0) * vv(2) * std::log(vv(1) + r);
                    tmpgrn += vv(0) * vv(1) * std::log(vv(2) + r);
                }
                else {
                    tmpgrn = -(vv(0) * vv(1) * vv(2) * (-9 * std::pow(interaction_radius_m, 2)
                             + std::pow(r, 2))) / (6 * std::pow(interaction_radius_m, 3));
                }

                tmpgreen_m[i][j][k] = tmpgrn / cellVolume;

            }
        }
    }

    
    Index I = nr_m[0] + 1;
    Index J = nr_m[1] + 1;
    Index K = nr_m[2] + 1;

    // the actual integration
    rho2_m = 0.0;
    rho2_m[I][J][K]  = tmpgreen_m[I+1][J+1][K+1];
    rho2_m[I][J][K] += tmpgreen_m[I+1][J][K];
    rho2_m[I][J][K] += tmpgreen_m[I][J+1][K];
    rho2_m[I][J][K] += tmpgreen_m[I][J][K+1];
    rho2_m[I][J][K] -= tmpgreen_m[I+1][J+1][K];
    rho2_m[I][J][K] -= tmpgreen_m[I+1][J][K+1];
    rho2_m[I][J][K] -= tmpgreen_m[I][J+1][K+1];
    rho2_m[I][J][K] -= tmpgreen_m[I][J][K];

    mirrorRhoField();

    fft_mp->transform(-1, rho2_m, grntr_m);

}

void P3MPoissonSolver::mirrorRhoField() {

    Index aI(0, 2 * nr_m[0]);
    Index aJ(0, 2 * nr_m[1]);

    Index J(0, nr_m[1]);
    Index K(0, nr_m[2]);

    Index IE(nr_m[0] + 1, 2 * nr_m[0]);
    Index JE(nr_m[1] + 1, 2 * nr_m[1]);
    Index KE(nr_m[2] + 1, 2 * nr_m[2]);

    Index mirroredIE = 2 * nr_m[0] - IE;
    Index mirroredJE = 2 * nr_m[1] - JE;
    Index mirroredKE = 2 * nr_m[2] - KE;

    rho2_m[0][0][0] = rho2_m[0][0][1];

    rho2_m[IE][J ][K ] = rho2_m[mirroredIE][J         ][K         ];
    rho2_m[aI][JE][K ] = rho2_m[aI        ][mirroredJE][K         ];
    rho2_m[aI][aJ][KE] = rho2_m[aI        ][aJ        ][mirroredKE];

}

Inform &P3MPoissonSolver::print(Inform& os) const {
    os << "* ************* P 3 M - P o i s s o n S o l v e r *************** " << endl;
    os << "* h        " << hr_m << '\n';
    os << "* RC       " << interaction_radius_m << '\n';
    os << "* ALPHA    " << alpha_m << '\n';
    os << "* *************************************************************** " << endl;
    return os;
}
