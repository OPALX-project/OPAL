//
// Class PartBunch
//   Particle Bunch.
//   A representation of a particle bunch as a vector of particles.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#include "Algorithms/PartBunch.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <cfloat>
#include <memory>
#include <utility>


#include "Distribution/Distribution.h"  // OPAL file
#include "Structure/FieldSolver.h"      // OPAL file
#include "Utilities/GeneralClassicException.h"

#include "Algorithms/ListElem.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_qrng.h>

#ifdef DBG_SCALARFIELD
    #include "Structure/FieldWriter.h"
#endif



//#define FIELDSTDOUT

PartBunch::PartBunch(const PartData *ref): // Layout is set using setSolver()
    PartBunchBase<double, 3>(new PartBunch::pbase_t(new Layout_t()), ref),
    interpolationCacheSet_m(false)
{

}


PartBunch::~PartBunch() {

}

// PartBunch::pbase_t* PartBunch::clone() {
//     return new pbase_t(new Layout_t());
// }


void PartBunch::initialize(FieldLayout_t *fLayout) {
    Layout_t* layout = static_cast<Layout_t*>(&getLayout());
    layout->getLayout().changeDomain(*fLayout);
}

void PartBunch::runTests() {

    Vector_t ll(-0.005);
    Vector_t ur(0.005);

    setBCAllPeriodic();

    NDIndex<3> domain = getFieldLayout().getDomain();
    for (unsigned int i = 0; i < Dimension; i++)
        nr_m[i] = domain[i].length();

    for (int i = 0; i < 3; i++)
        hr_m[i] = (ur[i] - ll[i]) / nr_m[i];

    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin(ll);

    rho_m.initialize(getMesh(),
                     getFieldLayout(),
                     GuardCellSizes<Dimension>(1),
                     bc_m);
    eg_m.initialize(getMesh(),
                    getFieldLayout(),
                    GuardCellSizes<Dimension>(1),
                    vbc_m);

    fs_m->solver_m->test(this);
}


void PartBunch::do_binaryRepart() {
    get_bounds(rmin_m, rmax_m);

    pbase_t* underlyingPbase =
        dynamic_cast<pbase_t*>(pbase.get());

    BinaryRepartition(*underlyingPbase);
    update();
    get_bounds(rmin_m, rmax_m);
    boundp();
}


void PartBunch::computeSelfFields(int binNumber) {
    IpplTimings::startTimer(selfFieldTimer_m);

    /// Set initial charge density to zero. Create image charge
    /// potential field.
    rho_m = 0.0;
    Field_t imagePotential = rho_m;

    /// Set initial E field to zero.
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {
        /// Mesh the whole domain
        if(fs_m->getFieldSolverType() == "SAAMG")
            resizeMesh();

        /// Scatter charge onto space charge grid.
        this->Q *= this->dt;
        if(!interpolationCacheSet_m) {
            if(interpolationCache_m.size() < getLocalNum()) {
                interpolationCache_m.create(getLocalNum() - interpolationCache_m.size());
            } else {
                interpolationCache_m.destroy(interpolationCache_m.size() - getLocalNum(),
                                             getLocalNum(),
                                             true);
            }
            interpolationCacheSet_m = true;
            this->Q.scatter(this->rho_m, this->R, IntrplCIC_t(), interpolationCache_m);
        } else {
            this->Q.scatter(this->rho_m, IntrplCIC_t(), interpolationCache_m);
        }

        this->Q /= this->dt;
        this->rho_m /= getdT();

        /// Calculate mesh-scale factor and get gamma for this specific slice (bin).
        double scaleFactor = 1;
        // double scaleFactor = Physics::c * getdT();
        double gammaz = getBinGamma(binNumber);

        /// Scale charge density to get charge density in real units. Account for
        /// Lorentz transformation in longitudinal direction.
        double tmp2 = 1 / hr_m[0] * 1 / hr_m[1] * 1 / hr_m[2] / (scaleFactor * scaleFactor * scaleFactor) / gammaz;
        rho_m *= tmp2;

        /// Scale mesh spacing to real units (meters). Lorentz transform the
        /// longitudinal direction.
        Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
        hr_scaled[2] *= gammaz;

        /// Find potential from charge in this bin (no image yet) using Poisson's
        /// equation (without coefficient: -1/(eps)).
        imagePotential = rho_m;

        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        /// Scale mesh back (to same units as particle locations.)
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// The scalar potential is given back in rho_m
        /// and must be converted to the right units.
        rho_m *= getCouplingConstant();

        /// IPPL Grad numerical computes gradient to find the
        /// electric field (in bin rest frame).
        eg_m = -Grad(rho_m, eg_m);

        /// Scale field. Combine Lorentz transform with conversion to proper field
        /// units.
        eg_m *= Vector_t(gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        // If desired write E-field and potential to terminal
#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        int mx = (int)nr_m[0];
        int mx2 = (int)nr_m[0] / 2;
        int my = (int)nr_m[1];
        int my2 = (int)nr_m[1] / 2;
        int mz = (int)nr_m[2];
        int mz2 = (int)nr_m[2] / 2;

        for (int i=0; i<mx; i++ )
            *gmsg << "Bin " << binNumber
                  << ", Self Field along x axis E = " << eg_m[i][my2][mz2]
                  << ", Pot = " << rho_m[i][my2][mz2]  << endl;

        for (int i=0; i<my; i++ )
            *gmsg << "Bin " << binNumber
                  << ", Self Field along y axis E = " << eg_m[mx2][i][mz2]
                  << ", Pot = " << rho_m[mx2][i][mz2]  << endl;

        for (int i=0; i<mz; i++ )
            *gmsg << "Bin " << binNumber
                  << ", Self Field along z axis E = " << eg_m[mx2][my2][i]
                  << ", Pot = " << rho_m[mx2][my2][i]  << endl;
#endif

        /// Interpolate electric field at particle positions.  We reuse the
        /// cached information about where the particles are relative to the
        /// field, since the particles have not moved since this the most recent
        /// scatter operation.
        Eftmp.gather(eg_m, IntrplCIC_t(), interpolationCache_m);
        //Eftmp.gather(eg_m, this->R, IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the electric field.
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
         *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
         *  \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        double betaC = std::sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;

        Bf(0) = Bf(0) - betaC * Eftmp(1);
        Bf(1) = Bf(1) + betaC * Eftmp(0);

        Ef += Eftmp;

        /// Now compute field due to image charge. This is done separately as the image charge
        /// is moving to -infinity (instead of +infinity) so the Lorentz transform is different.

        /// Find z shift for shifted Green's function.
        NDIndex<3> domain = getFieldLayout().getDomain();
        Vector_t origin = rho_m.get_mesh().get_origin();
        double hz = rho_m.get_mesh().get_meshSpacing(2);
        double zshift = -(2 * origin(2) + (domain[2].first() + domain[2].last() + 1) * hz) * gammaz * scaleFactor;

        /// Find potential from image charge in this bin using Poisson's
        /// equation (without coefficient: -1/(eps)).
        fs_m->solver_m->computePotential(imagePotential, hr_scaled, zshift);

        /// Scale mesh back (to same units as particle locations.)
        imagePotential *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// The scalar potential is given back in rho_m
        /// and must be converted to the right units.
        imagePotential *= getCouplingConstant();

#ifdef DBG_SCALARFIELD
        const int dumpFreq = 100;
        VField_t tmp_eg = eg_m;

        if ((localTrackStep_m + 1) % dumpFreq == 0) {
            FieldWriter fwriter;
            fwriter.dumpField(rho_m, "phi", "V", localTrackStep_m / dumpFreq, &imagePotential);
        }
#endif

        /// IPPL Grad numerical computes gradient to find the
        /// electric field (in rest frame of this bin's image
        /// charge).
        eg_m = -Grad(imagePotential, eg_m);

        /// Scale field. Combine Lorentz transform with conversion to proper field
        /// units.
        eg_m *= Vector_t(gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        // If desired write E-field and potential to terminal
#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        //int mx = (int)nr_m[0];
        //int mx2 = (int)nr_m[0] / 2;
        //int my = (int)nr_m[1];
        //int my2 = (int)nr_m[1] / 2;
        //int mz = (int)nr_m[2];
        //int mz2 = (int)nr_m[2] / 2;

        for (int i=0; i<mx; i++ )
            *gmsg << "Bin " << binNumber
                  << ", Image Field along x axis E = " << eg_m[i][my2][mz2]
                  << ", Pot = " << rho_m[i][my2][mz2]  << endl;

        for (int i=0; i<my; i++ )
            *gmsg << "Bin " << binNumber
                  << ", Image Field along y axis E = " << eg_m[mx2][i][mz2]
                  << ", Pot = " << rho_m[mx2][i][mz2]  << endl;

        for (int i=0; i<mz; i++ )
            *gmsg << "Bin " << binNumber
                  << ", Image Field along z axis E = " << eg_m[mx2][my2][i]
                  << ", Pot = " << rho_m[mx2][my2][i]  << endl;
#endif

#ifdef DBG_SCALARFIELD
        tmp_eg += eg_m;
        if ((localTrackStep_m + 1) % dumpFreq == 0) {
            FieldWriter fwriter;
            fwriter.dumpField(tmp_eg, "e", "V/m", localTrackStep_m / dumpFreq);
        }
#endif

        /// Interpolate electric field at particle positions.  We reuse the
        /// cached information about where the particles are relative to the
        /// field, since the particles have not moved since this the most recent
        /// scatter operation.
        Eftmp.gather(eg_m, IntrplCIC_t(), interpolationCache_m);
        //Eftmp.gather(eg_m, this->R, IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the image charge electric field. Note that beta will have
         *  the opposite sign from the bunch charge field, as the image charge is moving in the opposite direction.
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
         *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
         *  \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        Bf(0) = Bf(0) + betaC * Eftmp(1);
        Bf(1) = Bf(1) - betaC * Eftmp(0);

        Ef += Eftmp;

    }
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::resizeMesh() {
    double xmin = fs_m->solver_m->getXRangeMin();
    double xmax = fs_m->solver_m->getXRangeMax();
    double ymin = fs_m->solver_m->getYRangeMin();
    double ymax = fs_m->solver_m->getYRangeMax();
    double zmin = fs_m->solver_m->getZRangeMin();
    double zmax = fs_m->solver_m->getZRangeMax();

    if(xmin > rmin_m[0] || xmax < rmax_m[0] ||
       ymin > rmin_m[1] || ymax < rmax_m[1]) {

        for (unsigned int n = 0; n < getLocalNum(); n++) {

            if(R[n](0) < xmin || R[n](0) > xmax ||
               R[n](1) < ymin || R[n](1) > ymax) {

                // delete the particle
                INFOMSG(level2 << "destroyed particle with id=" << ID[n] << endl;);
                destroy(1, n);
            }

        }

        update();
        boundp();
        get_bounds(rmin_m, rmax_m);
    }
    Vector_t mymin = Vector_t(xmin, ymin , zmin);
    Vector_t mymax = Vector_t(xmax, ymax , zmax);

    for (int i = 0; i < 3; i++)
        hr_m[i]   = (mymax[i] - mymin[i])/nr_m[i];

    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin(mymin);

    rho_m.initialize(getMesh(),
                     getFieldLayout(),
                     GuardCellSizes<Dimension>(1),
                     bc_m);
    eg_m.initialize(getMesh(),
                    getFieldLayout(),
                    GuardCellSizes<Dimension>(1),
                    vbc_m);

    update();

//    setGridIsFixed();
}

void PartBunch::computeSelfFields() {
    IpplTimings::startTimer(selfFieldTimer_m);
    rho_m = 0.0;
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {
        //mesh the whole domain
        if(fs_m->getFieldSolverType() == "SAAMG")
            resizeMesh();

        //scatter charges onto grid
        this->Q *= this->dt;
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());
        this->Q /= this->dt;
        this->rho_m /= getdT();

        //calculating mesh-scale factor
        double gammaz = sum(this->P)[2] / getTotalNum();
        gammaz *= gammaz;
        gammaz = std::sqrt(gammaz + 1.0);
        double scaleFactor = 1;
        // double scaleFactor = Physics::c * getdT();
        //and get meshspacings in real units [m]
        Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
        hr_scaled[2] *= gammaz;

        //double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / (scaleFactor*scaleFactor*scaleFactor) / gammaz;
        double tmp2 = 1 / hr_scaled[0] * 1 / hr_scaled[1] * 1 / hr_scaled[2];
        //divide charge by a 'grid-cube' volume to get [C/m^3]
        rho_m *= tmp2;

#ifdef DBG_SCALARFIELD
        FieldWriter fwriter;
        fwriter.dumpField(rho_m, "rho", "C/m^3", localTrackStep_m);
#endif

        // charge density is in rho_m
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        //do the multiplication of the grid-cube volume coming
        //from the discretization of the convolution integral.
        //this is only necessary for the FFT solver
        //FIXME: later move this scaling into FFTPoissonSolver
        if(fs_m->getFieldSolverType() == "FFT" || fs_m->getFieldSolverType() == "FFTBOX")
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        // the scalar potential is given back in rho_m in units
        // [C/m] = [F*V/m] and must be divided by
        // 4*pi*\epsilon_0 [F/m] resulting in [V]
        rho_m *= getCouplingConstant();

        //write out rho
#ifdef DBG_SCALARFIELD
        fwriter.dumpField(rho_m, "phi", "V", localTrackStep_m);
#endif

        // IPPL Grad divides by hr_m [m] resulting in
        // [V/m] for the electric field
        eg_m = -Grad(rho_m, eg_m);

        eg_m *= Vector_t(gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        //write out e field
#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        int mx = (int)nr_m[0];
        int mx2 = (int)nr_m[0] / 2;
        int my = (int)nr_m[1];
        int my2 = (int)nr_m[1] / 2;
        int mz = (int)nr_m[2];
        int mz2 = (int)nr_m[2] / 2;

        for (int i=0; i<mx; i++ )
            *gmsg << "Field along x axis Ex = " << eg_m[i][my2][mz2] << " Pot = " << rho_m[i][my2][mz2]  << endl;

        for (int i=0; i<my; i++ )
            *gmsg << "Field along y axis Ey = " << eg_m[mx2][i][mz2] << " Pot = " << rho_m[mx2][i][mz2]  << endl;

        for (int i=0; i<mz; i++ )
            *gmsg << "Field along z axis Ez = " << eg_m[mx2][my2][i] << " Pot = " << rho_m[mx2][my2][i]  << endl;
#endif

#ifdef DBG_SCALARFIELD
        fwriter.dumpField(eg_m, "e", "V/m", localTrackStep_m);
#endif

        // interpolate electric field at particle positions.  We reuse the
        // cached information about where the particles are relative to the
        // field, since the particles have not moved since this the most recent
        // scatter operation.
        Ef.gather(eg_m, this->R,  IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the eletric field
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
         *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
         *  \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        double betaC = std::sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;

        Bf(0) = Bf(0) - betaC * Ef(1);
        Bf(1) = Bf(1) + betaC * Ef(0);
    }
    IpplTimings::stopTimer(selfFieldTimer_m);
}

/**
 * \method computeSelfFields_cycl()
 * \brief Calculates the self electric field from the charge density distribution for use in cyclotrons
 * \see ParallelCyclotronTracker
 * \warning none yet
 *
 * Comments -DW:
 * I have made some changes in here:
 * -) Some refacturing to make more similar to computeSelfFields()
 * -) Added meanR and quaternion to be handed to the function so that SAAMG solver knows how to rotate the boundary geometry correctly.
 * -) Fixed an error where gamma was not taken into account correctly in direction of movement (y in cyclotron)
 * -) Comment: There is no account for image charges in the cyclotron tracker (yet?)!
 */
void PartBunch::computeSelfFields_cycl(double gamma) {

    IpplTimings::startTimer(selfFieldTimer_m);

    /// set initial charge density to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {
        /// mesh the whole domain
        if(fs_m->getFieldSolverType() == "SAAMG")
            resizeMesh();

        /// scatter particles charge onto grid.
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length (y for cyclotron) enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// from charge (C) to charge density (C/m^3).
        double tmp2 = 1.0 / (hr_scaled[0] * hr_scaled[1] * hr_scaled[2]);
        rho_m *= tmp2;

        // If debug flag is set, dump scalar field (charge density 'rho') into file under ./data/
#ifdef DBG_SCALARFIELD
        FieldWriter fwriter;
        fwriter.dumpField(rho_m, "rho", "C/m^3", localTrackStep_m);
#endif

        /// now charge density is in rho_m
        /// calculate Possion equation (without coefficient: -1/(eps))
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        //do the multiplication of the grid-cube volume coming
        //from the discretization of the convolution integral.
        //this is only necessary for the FFT solver
        //TODO FIXME: later move this scaling into FFTPoissonSolver
        if(fs_m->getFieldSolverType() == "FFT" || fs_m->getFieldSolverType() == "FFTBOX")
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// retrive coefficient: -1/(eps)
        rho_m *= getCouplingConstant();

        // If debug flag is set, dump scalar field (potential 'phi') into file under ./data/
#ifdef DBG_SCALARFIELD
        fwriter.dumpField(rho_m, "phi", "V", localTrackStep_m);
#endif

        /// calculate electric field vectors from field potential
        eg_m = -Grad(rho_m, eg_m);

        /// Back Lorentz transformation
        /// CAVE: y coordinate needs 1/gamma factor because IPPL function Grad() divides by
        /// hr_m which is not scaled appropriately with Lorentz contraction in y direction
        /// only hr_scaled is! -DW
        eg_m *= Vector_t(gamma, 1.0 / gamma, gamma);

#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        int mx = (int)nr_m[0];
        int mx2 = (int)nr_m[0] / 2;
        int my = (int)nr_m[1];
        int my2 = (int)nr_m[1] / 2;
        int mz = (int)nr_m[2];
        int mz2 = (int)nr_m[2] / 2;

        for (int i=0; i<mx; i++ )
            *gmsg << "Field along x axis Ex = " << eg_m[i][my2][mz2] << " Pot = " << rho_m[i][my2][mz2]  << endl;

        for (int i=0; i<my; i++ )
            *gmsg << "Field along y axis Ey = " << eg_m[mx2][i][mz2] << " Pot = " << rho_m[mx2][i][mz2]  << endl;

        for (int i=0; i<mz; i++ )
            *gmsg << "Field along z axis Ez = " << eg_m[mx2][my2][i] << " Pot = " << rho_m[mx2][my2][i]  << endl;
#endif

#ifdef DBG_SCALARFIELD
        fwriter.dumpField(eg_m, "e", "V/m", localTrackStep_m);
#endif

        /// interpolate electric field at particle positions.
        Ef.gather(eg_m, this->R,  IntrplCIC_t());

        /// calculate coefficient
        // Relativistic E&M says gamma*v/c^2 = gamma*beta/c = sqrt(gamma*gamma-1)/c
        // but because we already transformed E_trans into the moving frame we have to
        // add 1/gamma so we are using the E_trans from the rest frame -DW
        double betaC = std::sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

        /// calculate B field from E field
        Bf(0) =  betaC * Ef(2);
        Bf(2) = -betaC * Ef(0);

    }

    /*
    *gmsg << "gamma =" << gamma << endl;
    *gmsg << "dx,dy,dz =(" << hr_m[0] << ", " << hr_m[1] << ", " << hr_m[2] << ") [m] " << endl;
    *gmsg << "max of bunch is (" << rmax_m(0) << ", " << rmax_m(1) << ", " << rmax_m(2) << ") [m] " << endl;
    *gmsg << "min of bunch is (" << rmin_m(0) << ", " << rmin_m(1) << ", " << rmin_m(2) << ") [m] " << endl;
    */

    IpplTimings::stopTimer(selfFieldTimer_m);
}

/**
 * \method computeSelfFields_cycl()
 * \brief Calculates the self electric field from the charge density distribution for use in cyclotrons
 * \see ParallelCyclotronTracker
 * \warning none yet
 *
 * Overloaded version for having multiple bins with separate gamma for each bin. This is necessary
 * For multi-bunch mode.
 *
 * Comments -DW:
 * I have made some changes in here:
 * -) Some refacturing to make more similar to computeSelfFields()
 * -) Added meanR and quaternion to be handed to the function (TODO: fall back to meanR = 0 and unit quaternion
 *    if not specified) so that SAAMG solver knows how to rotate the boundary geometry correctly.
 * -) Fixed an error where gamma was not taken into account correctly in direction of movement (y in cyclotron)
 * -) Comment: There is no account for image charges in the cyclotron tracker (yet?)!
 */
void PartBunch::computeSelfFields_cycl(int bin) {
    IpplTimings::startTimer(selfFieldTimer_m);

    /// set initial charge dentsity to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    /// get gamma of this bin
    double gamma = getBinGamma(bin);

    if(fs_m->hasValidSolver()) {
        /// mesh the whole domain
        if(fs_m->getFieldSolverType() == "SAAMG")
            resizeMesh();

        /// scatter particles charge onto grid.
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length (y for cyclotron) enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// from charge (C) to charge density (C/m^3).
        double tmp2 = 1.0 / (hr_scaled[0] * hr_scaled[1] * hr_scaled[2]);
        rho_m *= tmp2;

        // If debug flag is set, dump scalar field (charge density 'rho') into file under ./data/
#ifdef DBG_SCALARFIELD
        FieldWriter fwriter;
        fwriter.dumpField(rho_m, "rho", "C/m^3", localTrackStep_m);
#endif

        /// now charge density is in rho_m
        /// calculate Possion equation (without coefficient: -1/(eps))
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        // Do the multiplication of the grid-cube volume coming from the discretization of the convolution integral.
        // This is only necessary for the FFT solver. FIXME: later move this scaling into FFTPoissonSolver
        if(fs_m->getFieldSolverType() == "FFT" || fs_m->getFieldSolverType() == "FFTBOX")
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// retrive coefficient: -1/(eps)
        rho_m *= getCouplingConstant();

        // If debug flag is set, dump scalar field (potential 'phi') into file under ./data/
#ifdef DBG_SCALARFIELD
        fwriter.dumpField(rho_m, "phi", "V", localTrackStep_m);
#endif

        /// calculate electric field vectors from field potential
        eg_m = -Grad(rho_m, eg_m);

        /// Back Lorentz transformation
        /// CAVE: y coordinate needs 1/gamma factor because IPPL function Grad() divides by
        /// hr_m which is not scaled appropriately with Lorentz contraction in y direction
        /// only hr_scaled is! -DW
        eg_m *= Vector_t(gamma, 1.0 / gamma, gamma);

#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        int mx = (int)nr_m[0];
        int mx2 = (int)nr_m[0] / 2;
        int my = (int)nr_m[1];
        int my2 = (int)nr_m[1] / 2;
        int mz = (int)nr_m[2];
        int mz2 = (int)nr_m[2] / 2;

        for (int i=0; i<mx; i++ )
            *gmsg << "Bin " << bin
                  << ", Field along x axis Ex = " << eg_m[i][my2][mz2]
                  << ", Pot = " << rho_m[i][my2][mz2]  << endl;

        for (int i=0; i<my; i++ )
            *gmsg << "Bin " << bin
                  << ", Field along y axis Ey = " << eg_m[mx2][i][mz2]
                  << ", Pot = " << rho_m[mx2][i][mz2]  << endl;

        for (int i=0; i<mz; i++ )
            *gmsg << "Bin " << bin
                  << ", Field along z axis Ez = " << eg_m[mx2][my2][i]
                  << ", Pot = " << rho_m[mx2][my2][i]  << endl;
#endif

        // If debug flag is set, dump vector field (electric field) into file under ./data/
#ifdef DBG_SCALARFIELD
        fwriter.dumpField(eg_m, "e", "V/m", localTrackStep_m);
#endif

        /// Interpolate electric field at particle positions.
        Eftmp.gather(eg_m, this->R,  IntrplCIC_t());

        /// Calculate coefficient
        double betaC = std::sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

        /// Calculate B_bin field from E_bin field accumulate B and E field
        Bf(0) = Bf(0) + betaC * Eftmp(2);
        Bf(2) = Bf(2) - betaC * Eftmp(0);

        Ef += Eftmp;
    }

    /*
    *gmsg << "gamma =" << gamma << endl;
    *gmsg << "dx,dy,dz =(" << hr_m[0] << ", " << hr_m[1] << ", " << hr_m[2] << ") [m] " << endl;
    *gmsg << "max of bunch is (" << rmax_m(0) << ", " << rmax_m(1) << ", " << rmax_m(2) << ") [m] " << endl;
    *gmsg << "min of bunch is (" << rmin_m(0) << ", " << rmin_m(1) << ", " << rmin_m(2) << ") [m] " << endl;
    */


    IpplTimings::stopTimer(selfFieldTimer_m);
}


// void PartBunch::setMesh(Mesh_t* mesh) {
//     Layout_t* layout = static_cast<Layout_t*>(&getLayout());
// //     layout->getLayout().setMesh(mesh);
// }


// void PartBunch::setFieldLayout(FieldLayout_t* fLayout) {
//     Layout_t* layout = static_cast<Layout_t*>(&getLayout());
// //     layout->getLayout().setFieldLayout(fLayout);
// //     layout->rebuild_neighbor_data();
//     layout->getLayout().changeDomain(*fLayout);
// }


FieldLayout_t &PartBunch::getFieldLayout() {
    Layout_t* layout = static_cast<Layout_t*>(&getLayout());
    return dynamic_cast<FieldLayout_t &>(layout->getLayout().getFieldLayout());
}

void PartBunch::setBCAllPeriodic() {
    for (int i = 0; i < 2 * 3; ++i) {

        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dimension, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t, Dimension, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, Dimension, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t, Dimension, Mesh_t, Center_t>(i);
        }
        getBConds()[i] =  ParticlePeriodicBCond;
    }
    dcBeam_m=true;
    INFOMSG(level3 << "BC set P3M, all periodic" << endl);
}

void PartBunch::setBCAllOpen() {
    for (int i = 0; i < 2 * 3; ++i) {
        bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
        vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
        getBConds()[i] = ParticleNoBCond;
    }
    dcBeam_m=false;
    INFOMSG(level3 << "BC set for normal Beam" << endl);
}

void PartBunch::setBCForDCBeam() {
    for (int i = 0; i < 2 * 3; ++ i) {
        if (i >= 4) {
            if (Ippl::getNodes() > 1) {
                bc_m[i] = new ParallelPeriodicFace<double, 3, Mesh_t, Center_t>(i);
                vbc_m[i] = new ParallelPeriodicFace<Vector_t, 3, Mesh_t, Center_t>(i);
            } else {
                bc_m[i] = new PeriodicFace<double, 3, Mesh_t, Center_t>(i);
                vbc_m[i] = new PeriodicFace<Vector_t, 3, Mesh_t, Center_t>(i);
            }

            getBConds()[i] = ParticlePeriodicBCond;
        } else {
            bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
            getBConds()[i] = ParticleNoBCond;
        }
    }
    dcBeam_m=true;
    INFOMSG(level3 << "BC set for DC-Beam, longitudinal periodic" << endl);
}


void PartBunch::updateDomainLength(Vektor<int, 3>& grid) {
    NDIndex<3> domain = getFieldLayout().getDomain();
    for (unsigned int i = 0; i < Dimension; i++)
        grid[i] = domain[i].length();
}


void PartBunch::updateFields(const Vector_t& /*hr*/, const Vector_t& origin) {
    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin(origin);
    rho_m.initialize(getMesh(),
                     getFieldLayout(),
                     GuardCellSizes<Dimension>(1),
                     bc_m);
    eg_m.initialize(getMesh(),
                    getFieldLayout(),
                    GuardCellSizes<Dimension>(1),
                    vbc_m);
}


/**
 * Here we emit particles from the cathode. All particles in a new simulation (not a restart) initially reside in the bin
 container "pbin_m" and are not part of the beam bunch (so they cannot "see" fields, space charge etc.). In pbin_m, particles
 are sorted into the bins of a time histogram that describes the longitudinal time distribution of the beam, where the number
 of bins is given by \f$NBIN \times SBIN\f$. \f$NBIN\f$ and \f$SBIN\f$ are parameters given when defining the initial beam
 distribution. During emission, the time step of the simulation is set so that an integral number of these bins are emitted each step.
 Once all of the particles have been emitted, the simulation time step is reset to the value defined in the input file.

 A typical integration time step, \f$\Delta t\f$, is broken down into 3 sub-steps:

 1) Drift particles for \f$\frac{\Delta t}{2}\f$.

 2) Calculate fields and advance momentum.

 3) Drift particles for \f$\frac{\Delta t}{2}\f$ at the new momentum to complete the
 full time step.

 The difficulty for emission is that at the cathode position there is a step function discontinuity in the  fields. If we
 apply the typical integration time step across this boundary, we get an artificial numerical bunching of the beam, especially
 at very high accelerating fields. This function takes the cathode position boundary into account in order to achieve
 smoother particle emission.

 During an emission step, an integral number of time bins from the distribution histogram are emitted. However, each particle
 contained in those time bins will actually be emitted from the cathode at a different time, so will only spend some fraction
 of the total time step, \f$\Delta t_{full-timestep}\f$, in the simulation. The trick to emission is to give each particle
 a unique time step, \f$Delta t_{temp}\f$, that is equal to the actual time during the emission step that the particle
 exists in the simulation. For the next integration time step, the particle's time step is set back to the global time step,
 \f$\Delta t_{full-timestep}\f$.
  */


inline
PartBunch::VectorPair_t PartBunch::getEExtrema() {
    const Vector_t maxE = max(eg_m);
    //      const double maxL = max(dot(eg_m,eg_m));
    const Vector_t minE = min(eg_m);
    // INFOMSG("MaxE= " << maxE << " MinE= " << minE << endl);
    return VectorPair_t(maxE, minE);
}


inline
void PartBunch::resetInterpolationCache(bool clearCache) {
    interpolationCacheSet_m = false;
    if(clearCache) {
        interpolationCache_m.destroy(interpolationCache_m.size(), 0, true);
    }
}

void PartBunch::swap(unsigned int i, unsigned int j) {

    // FIXME
    PartBunchBase<double, 3>::swap(i, j);

    if (interpolationCacheSet_m)
        std::swap(interpolationCache_m[i], interpolationCache_m[j]);
}


Inform &PartBunch::print(Inform &os) {
    return PartBunchBase<double, 3>::print(os);
}