// -*- C++ -*-
/***************************************************************************
 *
 *
 * FFTBoxPoissonSolver.cc
 *
 *
 ***************************************************************************/

// include files
#include "Solvers/FFTBoxPoissonSolver.h"

#include "Algorithms/PartBunch.h"
#include "Physics/Physics.h"
#include "Utility/IpplTimings.h"
#include <fstream>
//////////////////////////////////////////////////////////////////////////////
// a little helper class to specialize the action of the Green's function
// calculation.  This should be specialized for each dimension
// to the proper action for computing the Green's function.  The first
// template parameter is the full type of the Field to compute, and the second
// is the dimension of the data, which should be specialized.

template<unsigned int Dim>
struct SpecializedGreensFunction { };

template<>
struct SpecializedGreensFunction<3> {
    template<class T, class FT, class FT2>
    static void calculate(Vektor<T, 3> &hrsq, FT &grn, FT2 *grnI) {
        grn = grnI[0] * hrsq[0] + grnI[1] * hrsq[1] + grnI[2] * hrsq[2];
        grn = 1.0 / sqrt(grn);
        grn[0][0][0] = grn[0][0][1];
    }
};

FFTBoxPoissonSolver::FFTBoxPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::string greensFunction, double boxSize):
    mesh_m(mesh),
    layout_m(fl),
    mesh2_m(0),
    layout2_m(0),
    greensFunction_m(greensFunction),
    a_m(boxSize)
{
    int i;
    domain_m = layout_m->getDomain();

    // For efficiency in the FFT's, we can use a parallel decomposition
    // which can be serial in the first dimension.
    // e_dim_tag decomp[3];
    // // e_dim_tag decomp2[3];
    // for(int d = 0; d < 3; ++d) {
    //     decomp[d] = layout_m->getRequestedDistribution(d);
    //     //decomp2[d] = layout_m->getRequestedDistribution(d);
    // }

    // The FFT's require double-sized field sizes in order to (more closely
    // do not understand this ...)
    // simulate an isolated system.  The FFT of the charge density field, rho,
    // would otherwise mimic periodic boundary conditions, i.e. as if there were
    // several beams set a periodic distance apart.  The double-sized fields
    // alleviate this problem.
    for(i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
        //domain2_m[i] = Index( 2*nr_m[i] + 1);
    }

    // create double sized mesh and layout objects for the use in the FFT's
    //mesh2_m = new Mesh_t(domain2_m);
    //layout2_m = new FieldLayout_t(*mesh2_m, decomp);
    //rho2_m.initialize(*mesh2_m, *layout2_m);
    grntr_m.initialize(*mesh_m, *layout_m);
    tmpgreen.initialize(*mesh_m, *layout_m);

    bool sineTransformDims[3];
    for(int d = 0; d < 3; ++d)
        sineTransformDims[d] = true;

    sine_m = new SINE_t(grntr_m.getDomain(), sineTransformDims, false);

    // these are fields that are used for calculating the Green's function.
    // they eliminate some calculation at each time-step.
    for(i = 0; i < 3; ++i) {
        grnIField_m[i].initialize(*mesh_m, *layout_m);
        grnIField_m[i][domain_m] = where(lt(domain_m[i], nr_m[i]),
                                         domain_m[i] * domain_m[i],
                                         (nr_m[i] - domain_m[i]) *
                                         (nr_m[i] - domain_m[i]));
    }

    GreensFunctionTimer_m = IpplTimings::getTimer("SF: GreenFTotal");

    if(greensFunction_m == std::string("INTEGRATED")) {
        IntGreensFunctionTimer1_m = IpplTimings::getTimer("SF: IntGreenF1");
        IntGreensFunctionTimer2_m = IpplTimings::getTimer("SF: IntGreenF2");
        IntGreensFunctionTimer3_m = IpplTimings::getTimer("SF: IntGreenF3");
        IntGreensFunctionTimer4_m = IpplTimings::getTimer("SF: IntGreenF4");

        ShIntGreensFunctionTimer1_m = IpplTimings::getTimer("SF: ShIntGreenF1");
        ShIntGreensFunctionTimer2_m = IpplTimings::getTimer("SF: ShIntGreenF2");
        ShIntGreensFunctionTimer3_m = IpplTimings::getTimer("SF: ShIntGreenF3");
        ShIntGreensFunctionTimer4_m = IpplTimings::getTimer("SF: ShIntGreenF4");
    } else {
        GreensFunctionTimer1_m = IpplTimings::getTimer("SF: GreenF1");
        GreensFunctionTimer4_m = IpplTimings::getTimer("SF: GreenF4");
    }
}

FFTBoxPoissonSolver::FFTBoxPoissonSolver(PartBunch &beam, std::string greensFunction):
    mesh_m(&beam.getMesh()),
    layout_m(&beam.getFieldLayout()),
    mesh2_m(0),
    layout2_m(0),
    greensFunction_m(greensFunction)
{
    int i;
    domain_m = layout_m->getDomain();

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
    for(i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
    }

    //rho_m.initialize(*mesh_m, *layout_m);
    grntr_m.initialize(*mesh_m, *layout_m);
    tmpgreen.initialize(*mesh_m, *layout_m);

    bool sineTransformDims[3];
    for(int d = 0; d < 3; ++d)
        sineTransformDims[d] = true;

    // create the FFT object
    sine_m = new SINE_t(grntr_m.getDomain(), sineTransformDims, false);

    // these are fields that are used for calculating the Green's function.
    // they eliminate some calculation at each time-step.
    for(i = 0; i < 3; ++i) {
        grnIField_m[i].initialize(*mesh_m, *layout_m);
        grnIField_m[i][domain_m] = where(lt(domain_m[i], nr_m[i]),
                                         domain_m[i] * domain_m[i],
                                         (nr_m[i] - domain_m[i]) *
                                         (nr_m[i] - domain_m[i]));
    }
    GreensFunctionTimer_m = IpplTimings::getTimer("SF: GreenFTotal");

    if(greensFunction_m == std::string("INTEGRATED")) {
        IntGreensFunctionTimer1_m = IpplTimings::getTimer("SF: IntGreenF1");
        IntGreensFunctionTimer2_m = IpplTimings::getTimer("SF: IntGreenF2");
        IntGreensFunctionTimer3_m = IpplTimings::getTimer("SF: IntGreenF3");
        IntGreensFunctionTimer4_m = IpplTimings::getTimer("SF: IntGreenF4");

        ShIntGreensFunctionTimer1_m = IpplTimings::getTimer("SF: ShIntGreenF1");
        ShIntGreensFunctionTimer2_m = IpplTimings::getTimer("SF: ShIntGreenF2");
        ShIntGreensFunctionTimer3_m = IpplTimings::getTimer("SF: ShIntGreenF3");
        ShIntGreensFunctionTimer4_m = IpplTimings::getTimer("SF: ShIntGreenF4");
    } else {
        GreensFunctionTimer1_m = IpplTimings::getTimer("SF: GreenF1");
        GreensFunctionTimer4_m = IpplTimings::getTimer("SF: GreenF4");
    }
}

////////////////////////////////////////////////////////////////////////////
// destructor
FFTBoxPoissonSolver::~FFTBoxPoissonSolver() {
    // delete the FFT object
    delete sine_m;
}

////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric potential from the image charge by solving
// the Poisson's equation

void FFTBoxPoissonSolver::computePotential(Field_t &rho, Vector_t hr, double zshift) {
    // use grid of complex doubled in both dimensions
    // and store rho in lower left quadrant of doubled grid
    //rho2_m = 0.0;

    //rho2_m[domain_m] = rho[domain_m];

    // needed in greens function
    hr_m = hr;
    // FFT double-sized charge density
    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    sine_m->transform(-1, rho);

    // must be called if the mesh size has changed
    // have to check if we can do G with h = (1,1,1)
    // and rescale later

    // Do image charge.
    // The minus sign is due to image charge.
    // Convolute transformed charge density with shifted green's function.
    IpplTimings::startTimer(GreensFunctionTimer_m);
    shiftedIntGreensFunction(zshift);
    IpplTimings::stopTimer(GreensFunctionTimer_m);

    // Multiply transformed charge density and
    // transformed Green's function.
    rho = - rho * grntr_m;

    // Inverse FFT to find image charge potential, rho2_m equals the electrostatic potential.
    sine_m->transform(+1, rho);

    // Re-use rho to store image potential. Flip z coordinate since this is a mirror image.
    Index I = nr_m[0];
    Index J = nr_m[1];
    Index K = nr_m[2];
    rho[I][J][K] = rho[I][J][nr_m[2] - K - 1];
}


////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation

void FFTBoxPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {
    // use grid of complex doubled in both dimensions
    // and store rho in lower left quadrant of doubled grid
    //rho2_m = 0.0;

    //rho2_m[domain_m] = rho[domain_m];

    // needed in greens function
    hr_m = hr;
    // FFT double-sized charge density
    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    sine_m->transform(-1, rho);
    // must be called if the mesh size has changed
    // have to check if we can do G with h = (1,1,1)
    // and rescale later
    IpplTimings::startTimer(GreensFunctionTimer_m);
    if(greensFunction_m == std::string("INTEGRATED"))
        integratedGreensFunction();
    else
        greensFunction();
    IpplTimings::stopTimer(GreensFunctionTimer_m);
    // multiply transformed charge density
    // and transformed Green function
    // Don't divide by (2*nx_m)*(2*ny_m), as Ryne does;
    // this normalization is done in POOMA's fft routine.
    rho *= grntr_m;

    // inverse FFT, rho2_m equals to the electrostatic potential
    sine_m->transform(+1, rho);
    // end convolution!

    // back to physical grid
    // reuse the charge density field to store the electrostatic potential
    //rho[domain_m] = rho2_m[domain_m];
}
///////////////////////////////////////////////////////////////////////////
// calculate the FFT of the Green's function for the given field
void FFTBoxPoissonSolver::greensFunction() {

    Vector_t hrsq(hr_m * hr_m);
    IpplTimings::startTimer(GreensFunctionTimer1_m);
    SpecializedGreensFunction<3>::calculate(hrsq, grntr_m, grnIField_m);
    IpplTimings::stopTimer(GreensFunctionTimer1_m);
    // Green's function calculation complete at this point.
    // The next step is to FFT it.

    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT

    IpplTimings::startTimer(GreensFunctionTimer4_m);
    sine_m->transform(-1, grntr_m);
    IpplTimings::stopTimer(GreensFunctionTimer4_m);
}

/** If the beam has a longitudinal size >> transverse size the
 * direct Green function at each mesh point is not efficient
 * (needs a lot of mesh points along the transverse size to
 * get a good resolution)
 *
 * If the charge density function is uniform within each cell
 * the following Green's function can be defined:
 *
 * \f[ \overline{G}(x_i - x_{i'}, y_j - y_{j'}, z_k - z_{k'}  cout << I << endl;
 cout << J << endl;
 cout << K << endl;
 cout << IE << endl;
 cout << JE << endl;
 cout << KE << endl;

 ) = \int_{x_{i'} - h_x/2}^{x_{i'} + h_x/2} dx' \int_{y_{j'} - h_y/2}^{y_{j'} + h_y/2} dy' \int_{z_{k'} - h_z/2}^{z_{k'} + h_z/2} dz' G(x_i - x_{i'}, y_j - y_{j'}, z_k - z_{k'}).
 * \f]
 */
void FFTBoxPoissonSolver::integratedGreensFunction() {

    tmpgreen = 0.0;
    double tmpgrn, r;
    NDIndex<3> idx =  layout2_m->getLocalNDIndex();

    IpplTimings::startTimer(IntGreensFunctionTimer1_m);

    /**
     * This integral can be calculated analytically in a closed from:
     */
    for(int k = idx[2].first(); k < std::min(nr_m[2] + 2, idx[2].last() + 3); k++) {
        for(int j = idx[1].first(); j < std::min(nr_m[1] + 2, idx[1].last() + 3); j++) {
            for(int i = idx[0].first(); i < std::min(nr_m[0] + 2, idx[0].last() + 3); i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = k * hr_m[2] - hr_m[2] / 2;

                r = sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                tmpgrn  = -vv(2) * vv(2) * atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                tmpgrn += -vv(1) * vv(1) * atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                tmpgrn += -vv(0) * vv(0) * atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                tmpgrn += vv(1) * vv(2) * log(vv(0) + r);
                tmpgrn += vv(0) * vv(2) * log(vv(1) + r);
                tmpgrn += vv(0) * vv(1) * log(vv(2) + r);

                tmpgreen[i][j][k] = tmpgrn / (hr_m[0] * hr_m[1] * hr_m[2]);

            }
        }
    }
    IpplTimings::stopTimer(IntGreensFunctionTimer1_m);

    IpplTimings::startTimer(IntGreensFunctionTimer2_m);

    grntr_m = 0.0;

    Index I = nr_m[0] + 1;
    Index J = nr_m[1] + 1;
    Index K = nr_m[2] + 1;

    // the actual integration
    grntr_m[I][J][K]  = tmpgreen[I+1][J+1][K+1];
    grntr_m[I][J][K] += tmpgreen[I+1][J][K];
    grntr_m[I][J][K] += tmpgreen[I][J+1][K];
    grntr_m[I][J][K] += tmpgreen[I][J][K+1];
    grntr_m[I][J][K] -= tmpgreen[I+1][J+1][K];
    grntr_m[I][J][K] -= tmpgreen[I+1][J][K+1];
    grntr_m[I][J][K] -= tmpgreen[I][J+1][K+1];
    grntr_m[I][J][K] -= tmpgreen[I][J][K];

    IpplTimings::stopTimer(IntGreensFunctionTimer2_m);
    IpplTimings::startTimer(IntGreensFunctionTimer3_m);

    //assign seems to have problems when we need values that are on another CPU, i.e. [I+1]
    /*assign(rho2_m[I][J][K] ,
      tmpgreen[I+1][J+1][K+1] - tmpgreen[I][J+1][K+1] - tmpgreen[I+1][J][K+1] + tmpgreen[I][J][K+1] - tmpgreen[I+1][J+1][K] +
      tmpgreen[I][J+1][K] + tmpgreen[I+1][J][K] - tmpgreen[I][J][K]);*/

    grntr_m[0][0][0] = grntr_m[0][0][1];

    Index IE(nr_m[0] + 1, 2 * nr_m[0] - 1);
    Index JE(nr_m[1] + 1, 2 * nr_m[1] - 1);
    Index KE(nr_m[2] + 1, 2 * nr_m[2] - 1);

    grntr_m[IE][JE][KE] = grntr_m[2*nr_m[0] - IE][2*nr_m[1] - JE][2*nr_m[2] - KE];
    grntr_m[IE][J ][K ] = grntr_m[2*nr_m[0] - IE][J][K];
    grntr_m[I ][JE][K ] = grntr_m[I][2*nr_m[1] - JE][K];
    grntr_m[I ][J ][KE] = grntr_m[I][J][2*nr_m[2] - KE];
    grntr_m[IE][JE][K ] = grntr_m[2*nr_m[0] - IE][2*nr_m[1] - JE][K];
    grntr_m[IE][J ][KE] = grntr_m[2*nr_m[0] - IE][J][2*nr_m[2] - KE];
    grntr_m[I ][JE][KE] = grntr_m[I][2*nr_m[1] - JE][2*nr_m[2] - KE];
    IpplTimings::stopTimer(IntGreensFunctionTimer3_m);

    IpplTimings::startTimer(IntGreensFunctionTimer4_m);
    sine_m->transform(-1, grntr_m);
    IpplTimings::stopTimer(IntGreensFunctionTimer4_m);
}


void FFTBoxPoissonSolver::shiftedIntGreensFunction(double zshift) {

    tmpgreen = 0.0;
    double tmpgrn, r;
    NDIndex<3> idx =  layout2_m->getLocalNDIndex();
    Field_t grn2, ggrn2;
    grn2.initialize(*mesh2_m, *layout2_m);
    grn2 = 0.0;
    ggrn2.initialize(*mesh2_m, *layout2_m);
    ggrn2 = 0.0;

    IpplTimings::startTimer(ShIntGreensFunctionTimer1_m);

    for(int k = idx[2].first(); k < std::min(nr_m[2] + 2, idx[2].last() + 3); k++) {
        for(int j = idx[1].first(); j < std::min(nr_m[1] + 2, idx[1].last() + 3); j++) {
            for(int i = idx[0].first(); i < std::min(nr_m[0] + 2, idx[0].last() + 3); i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = k * hr_m[2] - hr_m[2] / 2 + zshift;

                r = sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                tmpgrn  = -vv(2) * vv(2) * atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                tmpgrn += -vv(1) * vv(1) * atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                tmpgrn += -vv(0) * vv(0) * atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                tmpgrn += vv(1) * vv(2) * log(vv(0) + r);
                tmpgrn += vv(0) * vv(2) * log(vv(1) + r);
                tmpgrn += vv(0) * vv(1) * log(vv(2) + r);

                tmpgreen[i][j][k] = tmpgrn / (hr_m[0] * hr_m[1] * hr_m[2]);

            }
        }
    }
    IpplTimings::stopTimer(ShIntGreensFunctionTimer1_m);

    IpplTimings::startTimer(ShIntGreensFunctionTimer2_m);

    for(int k = idx[2].first(); k < std::min(nr_m[2] + 2, idx[2].last() + 3); k++) {
        for(int j = idx[1].first(); j < std::min(nr_m[1] + 2, idx[1].last() + 3); j++) {
            for(int i = idx[0].first(); i < std::min(nr_m[0] + 2, idx[0].last() + 3); i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = zshift - hr_m[2] * (nr_m[2] - k) - hr_m[2] / 2;
                //vv(2) = zshift-k*hr_m[2]-hr_m[2]/2;

                r = sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                tmpgrn  = -vv(2) * vv(2) * atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                tmpgrn += -vv(1) * vv(1) * atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                tmpgrn += -vv(0) * vv(0) * atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                tmpgrn += vv(1) * vv(2) * log(vv(0) + r);
                tmpgrn += vv(0) * vv(2) * log(vv(1) + r);
                tmpgrn += vv(0) * vv(1) * log(vv(2) + r);

                grn2[i][j][k] = tmpgrn / (hr_m[0] * hr_m[1] * hr_m[2]);

            }
        }
    }
    IpplTimings::stopTimer(ShIntGreensFunctionTimer2_m);

    IpplTimings::startTimer(ShIntGreensFunctionTimer3_m);

    grntr_m = 0.0;

    Index I = nr_m[0] + 1;
    Index J = nr_m[1] + 1;
    Index K = nr_m[2] + 1;

    // the actual integration
    grntr_m[I][J][K]  = tmpgreen[I+1][J+1][K+1];
    grntr_m[I][J][K] += tmpgreen[I+1][J][K];
    grntr_m[I][J][K] += tmpgreen[I][J+1][K];
    grntr_m[I][J][K] += tmpgreen[I][J][K+1];
    grntr_m[I][J][K] -= tmpgreen[I+1][J+1][K];
    grntr_m[I][J][K] -= tmpgreen[I+1][J][K+1];
    grntr_m[I][J][K] -= tmpgreen[I][J+1][K+1];
    grntr_m[I][J][K] -= tmpgreen[I][J][K];
    tmpgreen = 0.0;
    //later replace ggrn2 with tmpgreen
    ggrn2[I][J][K]  = grn2[I+1][J+1][K+1];
    ggrn2[I][J][K] += grn2[I+1][J][K];
    ggrn2[I][J][K] += grn2[I][J+1][K];
    ggrn2[I][J][K] += grn2[I][J][K+1];
    ggrn2[I][J][K] -= grn2[I+1][J+1][K];
    ggrn2[I][J][K] -= grn2[I+1][J][K+1];
    ggrn2[I][J][K] -= grn2[I][J+1][K+1];
    ggrn2[I][J][K] -= grn2[I][J][K];

    Index IE(nr_m[0] + 1, 2 * nr_m[0] - 1);
    Index JE(nr_m[1] + 1, 2 * nr_m[1] - 1);
    Index KE(nr_m[2], 2 * nr_m[2] - 1);

    /**
     ** (x[0:nr_m[0]-1]^2 + y[0:nr_m[1]-1]^2 + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[0:nr_m[1]-1]^2 + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     ** (x[0:nr_m[0]-1]^2 + y[nr_m[1]:1]^2   + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[nr_m[1]:1]^2   + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     **
     ** (x[0:nr_m[0]-1]^2 + y[0:nr_m[1]-1]^2 + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[0:nr_m[1]-1]^2 + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     ** (x[0:nr_m[0]-1]^2 + y[nr_m[1]:1]^2   + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[nr_m[1]:1]^2   + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     */

    grntr_m[IE][J ][K ] = grntr_m[2*nr_m[0] - IE][J][K];
    grntr_m[I ][JE][K ] = grntr_m[I][2*nr_m[1] - JE][K];
    grntr_m[IE][JE][K ] = grntr_m[2*nr_m[0] - IE][2*nr_m[1] - JE][K];

    grntr_m[I ][J ][KE] = ggrn2[I][J][KE - nr_m[2]];
    grntr_m[IE][J ][KE] = ggrn2[2*nr_m[0] - IE][J][KE - nr_m[2]];
    grntr_m[I ][JE][KE] = ggrn2[I][2*nr_m[1] - JE][KE - nr_m[2]];
    grntr_m[IE][JE][KE] = ggrn2[2*nr_m[0] - IE][2*nr_m[1] - JE][KE - nr_m[2]];

    IpplTimings::stopTimer(ShIntGreensFunctionTimer3_m);

    IpplTimings::startTimer(ShIntGreensFunctionTimer4_m);
    sine_m->transform(-1, grntr_m);
    IpplTimings::stopTimer(ShIntGreensFunctionTimer4_m);

}

Inform &FFTBoxPoissonSolver::print(Inform &os) const {
    os << "* ************* F F T P o i s s o n S o l v e r ************************************ " << endl;
    os << "* h " << hr_m << '\n';
    os << "* ********************************************************************************** " << endl;
    return os;
}

/***************************************************************************
 * $RCSfile: FFTBoxPoissonSolver.cc,v $   $Author: adelmann $
 * $Revision: 1.6 $   $Date: 2001/08/16 09:36:08 $
 ***************************************************************************/