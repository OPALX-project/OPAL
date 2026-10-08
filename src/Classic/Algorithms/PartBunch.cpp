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

#include <algorithm>
#include <cfloat>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Particle/ParticleBalancer.h"

#include "Algorithms/ListElem.h"
#include "Distribution/Distribution.h"
#include "Structure/FieldSolver.h"
#include "Utilities/GeneralClassicException.h"

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


void PartBunch::do_binaryRepart() {
    get_bounds(rmin_m, rmax_m);

    pbase_t* underlyingPbase =
        dynamic_cast<pbase_t*>(pbase_m.get());

    BinaryRepartition(*underlyingPbase);
    update();
    get_bounds(rmin_m, rmax_m);
    boundp();
}


void PartBunch::computeSelfFields(int binNumber) {
    IpplTimings::startTimer(selfFieldTimer_m);

    if (fs_m->getFieldSolverType() == FieldSolverType::P3M) {
        throw GeneralClassicException("PartBunch::computeSelfFields(int binNumber)", 
                            "P3M solver not available during emission");
    }

    /// Set initial charge density to zero. Create image charge
    /// potential field.
    rho_m = 0.0;
    Field_t imagePotential = rho_m;

    /// Set initial E field to zero.
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {
        /// Mesh the whole domain
        resizeMesh();

        /// Scatter charge onto space charge grid.
        checkUnmaskedInterpolation("PartBunch::computeSelfFields(int binNumber)");
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
        eg_m *= Vector_t({gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz)});

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
        eg_m *= Vector_t({gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz)});

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
    if (fs_m->getFieldSolverType() != FieldSolverType::SAAMG) {
        return;
    }

    double xmin = fs_m->solver_m->getXRangeMin();
    double xmax = fs_m->solver_m->getXRangeMax();
    double ymin = fs_m->solver_m->getYRangeMin();
    double ymax = fs_m->solver_m->getYRangeMax();

    if(xmin > rmin_m[0] || xmax < rmax_m[0] ||
       ymin > rmin_m[1] || ymax < rmax_m[1]) {

        for (unsigned int n = 0; n < getLocalNum(); n++) {

            if(R[n](0) < xmin || R[n](0) > xmax ||
               R[n](1) < ymin || R[n](1) > ymax) {

                // delete the particle
                INFOMSG(level2 << "destroyed particle with id=" << ID[n] << endl);
                destroy(1, n);
            }

        }

        update();
        boundp();
        get_bounds(rmin_m, rmax_m);
    }

    Vector_t origin = Vector_t({0.0, 0.0, 0.0});

    // update the mesh origin and mesh spacing hr_m
    fs_m->solver_m->resizeMesh(origin, hr_m, rmin_m, rmax_m, dh_m);

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

    update();

//    setGridIsFixed();
}

void PartBunch::computeSelfFields() {
    IpplTimings::startTimer(selfFieldTimer_m);
    rho_m = 0.0;
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {
        //mesh the whole domain
        resizeMesh();

        //scatter charges onto grid
        checkUnmaskedInterpolation("PartBunch::computeSelfFields()");
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

        double Npoints = nr_m[0] * nr_m[1] * nr_m[2];
        rmsDensity_m = std::sqrt((1.0 /Npoints) * sum((rho_m / Physics::q_e) * (rho_m / Physics::q_e)));

        calcDebyeLength(); 


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
        if (fs_m->getFieldSolverType() == FieldSolverType::FFT ||
            fs_m->getFieldSolverType() == FieldSolverType::FFTBOX) {
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];
        }

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

        if(fs_m->getFieldSolverType() == FieldSolverType::P3M) {
            fs_m->solver_m->calculatePairForces(this,gammaz);
        }

        Ef = Ef * Vector_t({gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz)});
        
        /** Magnetic field in x and y direction induced by the electric field
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

    if (getTotalNum() == 0) {
        return;
    }

    if (fs_m->getFieldSolverType() == FieldSolverType::SAAMG) {
        std::size_t numGridPoints = static_cast<std::size_t>(nr_m[0]) * static_cast<std::size_t>(nr_m[1]) * static_cast<std::size_t>(nr_m[2]);
        if (getTotalNum() < numGridPoints) {
            WARNMSG("The number of particles decreased below the grid points. Skip space-charge calculation." << endl);
            return;
        }
    }

    IpplTimings::startTimer(selfFieldTimer_m);

    if (fs_m->getFieldSolverType() == FieldSolverType::P3M) {
        throw GeneralClassicException("PartBunch::computeSelfFields_cycl(double gamma)",
                                      "P3M solver not available yet for cyclotrons");
    }

    /// set initial charge density to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    if (fs_m->hasValidSolver()) {
        /// mesh the whole domain
        resizeMesh();

        /// FIELDSOLVER, MESHFIT=CORE: the mesh of this solve, as the stat file and the log
        /// report it. A repartition after the fit of the solve step fits the mesh again.
        if (meshFit_m.type == MeshFitType::CORE) {
            scStatistics_m.meshSpacing = hr_m;
        }

        /// FIELDSOLVER, MESHFIT=CORE: on a mesh fitted to the core only the core particles
        /// are deposited, and only the particles within the range of cell centres gathered
        /// (CoreFitSC::ParticleClass); the others may lie outside the mesh.
        const bool coreMesh = isMeshFitToCore();
        std::vector<CoreFitSC::ParticleClass> particleClass;
        ParticleMask_t isCore;
        if (coreMesh) {
            const size_t localNum = getLocalNum();
            const ParticleAttrib<Vector_t>& Rc = this->R;
            const ParticleAttrib<int>& Binc = this->Bin;
            particleClass.resize(localNum);
            for (size_t i = 0; i < localNum; ++i) {
                particleClass[i] = CoreFitSC::classify(Binc[i], Rc[i],
                                                       scBoxMin_m, scBoxMax_m,
                                                       scMeshMin_m, scMeshMax_m);
            }
            isCore = [&particleClass](size_t i) {
                return particleClass[i] == CoreFitSC::ParticleClass::CORE;
            };
        } else {
            checkUnmaskedInterpolation("PartBunch::computeSelfFields_cycl(double gamma)");
        }

        /// scatter particles charge onto grid, on a mesh fitted to the core that of the core
        /// particles. One call serves both meshes (see scatterMasked()).
        scatterMasked(isCore);

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length (y for cyclotron) enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// from charge (C) to charge density (C/m^3).
        double tmp2 = 1.0 / (hr_scaled[0] * hr_scaled[1] * hr_scaled[2]);
        rho_m *= tmp2;

        // Diagnostics for the stat file only: a field reduction and two particle passes
        // with collectives, skipped when the tracker says no stat row needs them. On a mesh
        // fitted to the core (MESHFIT=CORE) the rms density, and with it the Debye length,
        // are those of the core charge on the core mesh.
        if (computeSCDiagnostics_m) {
            double Npoints = nr_m[0] * nr_m[1] * nr_m[2];
            rmsDensity_m = std::sqrt((1.0 /Npoints) * sum((rho_m / Physics::q_e) * (rho_m / Physics::q_e)));

            calcDebyeLength();
        }

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
        if (fs_m->getFieldSolverType() == FieldSolverType::FFT ||
            fs_m->getFieldSolverType() == FieldSolverType::FFTBOX) {
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];
        }

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
        eg_m *= Vector_t({gamma, 1.0 / gamma, gamma});

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

        if (coreMesh) {
            /// interpolate electric field at the particles within the range of cell
            /// centres, then the field of the particles outside the core bounds.
            gatherMasked([&particleClass](size_t i) {
                return particleClass[i] != CoreFitSC::ParticleClass::FAR;
            });
            computeFarField(gamma, particleClass);
        } else {
            /// interpolate electric field at particle positions.
            Ef.gather(eg_m, this->R,  IntrplCIC_t());
        }

        /// Calculate coefficient
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

    if (fs_m->getFieldSolverType() == FieldSolverType::P3M) {
        throw GeneralClassicException("PartBunch::computeSelfFields_cycl(int bin)", 
                            "P3M solver not available yet for cyclotrons");
    }

    /// set initial charge dentsity to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    /// get gamma of this bin
    double gamma = getBinGamma(bin);

    if(fs_m->hasValidSolver()) {
        /// mesh the whole domain
        resizeMesh();

        /// scatter particles charge onto grid.
        checkUnmaskedInterpolation("PartBunch::computeSelfFields_cycl(int bin)");
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
        if (fs_m->getFieldSolverType() == FieldSolverType::FFT ||
            fs_m->getFieldSolverType() == FieldSolverType::FFTBOX) {
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];
        }

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
        eg_m *= Vector_t({gamma, 1.0 / gamma, gamma});

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


void PartBunch::scatterMasked(const ParticleMask_t& accept) {
    /// Particles flagged lost (Bin < 0) stay in the bunch until the tracker deletes them,
    /// every DELPARTFREQ steps. Leave their charge off the mesh, as setBinCharge() does in
    /// the multi-bunch solve. Without a predicate their Q is set to zero for the scatter and
    /// restored after it: the integrator, ScatteringPhysics and the loss records still read
    /// it. With a predicate the accepted particles that are not flagged lost are copied, in
    /// order, into maskedQ_m and maskedR_m; a flagged particle adds only zero to the scatter
    /// of Q, so an all-true predicate gives the same rho_m, bit for bit. Both cases share
    /// one call of IPPL's scatter: with a second call, or with a loop of its own over IPPL's
    /// CIC kernel, the production build no longer inlines the kernel into
    /// computeSelfFields_cycl(double) and computeSelfFields().
    std::vector<std::pair<size_t, double>> lostCharge;
    if (accept) {
        const size_t localNum = getLocalNum();
        if (maskedQ_m.size() < localNum) {
            const size_t missing = localNum - maskedQ_m.size();
            maskedQ_m.create(missing);
            maskedR_m.create(missing);
        }
        size_t numAccepted = 0;
        for (size_t i = 0; i < localNum; ++i) {
            if (this->Bin[i] >= 0 && accept(i)) {
                maskedQ_m[numAccepted] = this->Q[i];
                maskedR_m[numAccepted] = this->R[i];
                ++numAccepted;
            }
        }
        const size_t numUnused = maskedQ_m.size() - numAccepted;
        maskedQ_m.destroy(numUnused, numAccepted);
        maskedR_m.destroy(numUnused, numAccepted);
    } else {
        for (size_t i = 0; i < getLocalNum(); ++i) {
            if (this->Bin[i] < 0) {
                lostCharge.emplace_back(i, this->Q[i]);
                this->Q[i] = 0.0;
            }
        }
    }

    const ParticleAttrib<double>& charge = accept ? maskedQ_m : this->Q;
    const ParticleAttrib<Vector_t>& position = accept ? maskedR_m : this->R;
    charge.scatter(rho_m, position, IntrplCIC_t());

    for (const auto& lost : lostCharge) {
        this->Q[lost.first] = lost.second;
    }
}


void PartBunch::gatherMasked(const ParticleMask_t& accept) {
    /// The loop of ParticleAttrib::gather, restricted to the accepted particles.
    eg_m.Uncompress();
    if (eg_m.isDirty()) {
        eg_m.fillGuardCells(true);
    }

    const Mesh_t& mesh = eg_m.get_mesh();
    for (size_t i = 0; i < getLocalNum(); ++i) {
        if (accept(i)) {
            IntrplCIC_t::gather(Ef[i], eg_m, this->R[i], mesh);
        }
    }

    eg_m.Compress();
}


void PartBunch::computeFarField(double gamma,
                                const std::vector<CoreFitSC::ParticleClass>& particleClass) {
    IpplTimings::startTimer(farFieldTimer_m);

    typedef CoreFitSC::ParticleClass Class_t;
    const int numNodes = Ippl::getNodes();
    const int myNode = Ippl::myNode();
    MPI_Comm comm = Ippl::getComm();
    const size_t localNum = getLocalNum();
    const ParticleAttrib<Vector_t>& Rc = this->R;
    const ParticleAttrib<double>& Qc = this->Q;
    const ParticleAttrib<int>& Binc = this->Bin;
    const ParticleIndex_t& IDc = this->ID;

    // rest-frame positions about the centre of the core bounds, for conditioning
    const Vector_t centre = 0.5 * (scBoxMin_m + scBoxMax_m);

    // This rank's part of the far table: the particles outside the core bounds, with their
    // rest-frame position, charge (0 for those flagged lost), class (1 for FAR, 0 for MESH)
    // and ID, which picks the subsample of the far-far sources.
    const int tableWidth = 6;
    std::vector<size_t> farIndex;
    std::vector<double> localTable;
    for (size_t i = 0; i < localNum; ++i) {
        if (particleClass[i] != Class_t::MESH && particleClass[i] != Class_t::FAR) {
            continue;
        }
        const Vector_t s = CoreFitSC::toRestFrame(Rc[i], centre, gamma);
        farIndex.push_back(i);
        localTable.insert(localTable.end(),
                          {s(0), s(1), s(2), Binc[i] >= 0 ? Qc[i] : 0.0,
                           particleClass[i] == Class_t::FAR ? 1.0 : 0.0,
                           static_cast<double>(IDc[i])});
    }

    // the table of all ranks, in rank order
    const int localCount = farIndex.size();
    std::vector<int> counts(numNodes), offsets(numNodes + 1, 0);
    MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
    for (int p = 0; p < numNodes; ++p) {
        offsets[p + 1] = offsets[p] + counts[p];
    }
    const size_t numTable = offsets[numNodes];
    scStatistics_m.numFlagged = 0;
    if (numTable == 0) {
        IpplTimings::stopTimer(farFieldTimer_m);
        return;
    }

    // monopole, centroid and quadrupole of the charge on the mesh
    CoreFitSC::MomentSums sums(centre, gamma);
    for (size_t i = 0; i < localNum; ++i) {
        if (particleClass[i] == Class_t::CORE) {
            sums.add(Qc[i], Rc[i]);
        }
    }
    const CoreFitSC::CoreMoments moments = sums.reduce(CoreFitSC::IpplReducer());

    std::vector<int> recvCounts(numNodes), displs(numNodes);
    for (int p = 0; p < numNodes; ++p) {
        recvCounts[p] = tableWidth * counts[p];
        displs[p] = tableWidth * offsets[p];
    }
    std::vector<double> table(tableWidth * numTable);
    MPI_Allgatherv(localTable.data(), tableWidth * localCount, MPI_DOUBLE,
                   table.data(), recvCounts.data(), displs.data(), MPI_DOUBLE, comm);

    // The live entries are the sources and the targets of the far-far field. The entries of
    // particles flagged lost get the multipole field (FAR) or keep their mesh field (MESH).
    std::vector<Vector_t> position(numTable);
    std::vector<size_t> liveRow;
    std::vector<Vector_t> livePosition;
    for (size_t j = 0; j < numTable; ++j) {
        const double* entry = &table[tableWidth * j];
        position[j] = Vector_t({entry[0], entry[1], entry[2]});
        if (entry[3] != 0.0) {
            liveRow.push_back(j);
            livePosition.push_back(position[j]);
        }
    }
    const size_t numLive = liveRow.size();

    // This rank's share of the table and of its live entries, by index: the far particles
    // sit on the few ranks that own the faces of the mesh nearest to them.
    const std::pair<size_t, size_t> range = CoreFitSC::getTableRange(numTable, myNode, numNodes);
    const std::pair<size_t, size_t> liveRange = CoreFitSC::getTableRange(numLive, myNode,
                                                                         numNodes);
    const FarFieldModel model = meshFit_m.farField;
    std::vector<Vector_t> multipole(numTable);
    std::vector<double> indicator(numTable);
    CoreFitSC::computeMultipoleField(moments, position, range.first, range.second,
                                     model != FarFieldModel::MONOPOLE, multipole, indicator);
    std::vector<Vector_t> liveFarFar(numLive, Vector_t(0.0));
    size_t farFarStride = 1;
    if (model == FarFieldModel::FULL) {
        // Above CoreFitSC::maxFarFarPairs pairs the sources are the live entries whose ID is
        // in the subsample, with their charge scaled to that of all live entries.
        farFarStride = CoreFitSC::getSampleStride(numLive, numLive, CoreFitSC::maxFarFarPairs);
        std::vector<Vector_t> sourcePosition;
        std::vector<double> sourceCharge;
        double liveCharge = 0.0;
        double sampledCharge = 0.0;
        for (const size_t j: liveRow) {
            const double* entry = &table[tableWidth * j];
            liveCharge += entry[3];
            if (CoreFitSC::isSampled(static_cast<size_t>(entry[5]), farFarStride)) {
                sourcePosition.push_back(position[j]);
                sourceCharge.push_back(entry[3]);
                sampledCharge += entry[3];
            }
        }
        if (farFarStride > 1 && sampledCharge != 0.0) {
            const double scale = liveCharge / sampledCharge;
            for (double& q: sourceCharge) {
                q *= scale;
            }
        }
        // softened by the smallest rest-frame cell of the core mesh
        const double eps = std::min({hr_m[0], gamma * hr_m[1], hr_m[2]});
        CoreFitSC::addDirectField(livePosition, liveRange.first, liveRange.second,
                                  sourcePosition, sourceCharge, eps, liveFarFar);
    }

    // The results of all ranks: the multipole field and the indicator of their share of the
    // table, then the far-far field of their share of the live entries.
    std::vector<double> localResult;
    localResult.reserve(4 * (range.second - range.first) +
                        3 * (liveRange.second - liveRange.first));
    for (size_t j = range.first; j < range.second; ++j) {
        localResult.insert(localResult.end(),
                           {multipole[j](0), multipole[j](1), multipole[j](2), indicator[j]});
    }
    for (size_t k = liveRange.first; k < liveRange.second; ++k) {
        localResult.insert(localResult.end(),
                           {liveFarFar[k](0), liveFarFar[k](1), liveFarFar[k](2)});
    }
    int numResults = 0;
    for (int p = 0; p < numNodes; ++p) {
        const std::pair<size_t, size_t> rangeP = CoreFitSC::getTableRange(numTable, p, numNodes);
        const std::pair<size_t, size_t> liveRangeP = CoreFitSC::getTableRange(numLive, p,
                                                                              numNodes);
        recvCounts[p] = 4 * (rangeP.second - rangeP.first) +
                        3 * (liveRangeP.second - liveRangeP.first);
        displs[p] = numResults;
        numResults += recvCounts[p];
    }
    std::vector<double> result(numResults);
    MPI_Allgatherv(localResult.data(), static_cast<int>(localResult.size()), MPI_DOUBLE,
                   result.data(), recvCounts.data(), displs.data(), MPI_DOUBLE, comm);
    std::vector<Vector_t> farFar(numTable, Vector_t(0.0));
    for (int p = 0; p < numNodes; ++p) {
        const std::pair<size_t, size_t> rangeP = CoreFitSC::getTableRange(numTable, p, numNodes);
        const std::pair<size_t, size_t> liveRangeP = CoreFitSC::getTableRange(numLive, p,
                                                                              numNodes);
        const double* entry = &result[displs[p]];
        for (size_t j = rangeP.first; j < rangeP.second; ++j, entry += 4) {
            multipole[j] = Vector_t({entry[0], entry[1], entry[2]});
            indicator[j] = entry[3];
        }
        for (size_t k = liveRangeP.first; k < liveRangeP.second; ++k, entry += 3) {
            farFar[liveRow[k]] = Vector_t({entry[0], entry[1], entry[2]});
        }
    }

    // The live FAR entries where the expansion is not valid get the exact field of the core
    // particles instead: every rank sums over its own, one sum over the ranks. Above
    // CoreFitSC::maxCorePairsPerCore pairs per core particle, the core particles in the outer
    // shell of the box and the others are summed apart (one more sum, of the shell count),
    // each over the particles whose ID is in its subsample if its pairs exceed that budget,
    // with their charge scaled to that of the part (one more sum).
    std::vector<size_t> flagged;
    for (const size_t j: liveRow) {
        const Class_t c = (table[tableWidth * j + 4] != 0.0 ? Class_t::FAR : Class_t::MESH);
        if (CoreFitSC::needsExactCoreField(c, indicator[j])) {
            flagged.push_back(j);
        }
    }
    scStatistics_m.numFlagged = flagged.size();
    // strides of the core particles in the shell of the box and of the others
    size_t coreStride[2] = {1, 1};
    if (model == FarFieldModel::FULL && !flagged.empty()) {
        IpplTimings::startTimer(exactCoreTimer_m);
        std::vector<Vector_t> target(flagged.size());
        for (size_t k = 0; k < flagged.size(); ++k) {
            target[k] = position[flagged[k]];
        }

        const double maxPairs = CoreFitSC::maxCorePairsPerCore * scNumCore_m;
        const bool split = (static_cast<double>(flagged.size()) * scNumCore_m > maxPairs);
        if (split) {
            double numShell = 0.0;
            for (size_t i = 0; i < localNum; ++i) {
                if (particleClass[i] == Class_t::CORE &&
                    CoreFitSC::isInCoreShell(Rc[i], scBoxMin_m, scBoxMax_m)) {
                    numShell += 1.0;
                }
            }
            allreduce(&numShell, 1, std::plus<double>());
            const size_t numShellCore = static_cast<size_t>(numShell);
            coreStride[0] = CoreFitSC::getSampleStride(flagged.size(), numShellCore, maxPairs);
            coreStride[1] = CoreFitSC::getSampleStride(flagged.size(),
                                                       scNumCore_m - numShellCore, maxPairs);
        }

        // the sources of the two parts (all of them in the first without a split), and the
        // charge of each part and of its subsample
        std::vector<Vector_t> corePosition[2];
        std::vector<double> coreCharge[2];
        double charges[4] = {0.0, 0.0, 0.0, 0.0};
        for (size_t i = 0; i < localNum; ++i) {
            if (particleClass[i] != Class_t::CORE) {
                continue;
            }
            const int part = (split && !CoreFitSC::isInCoreShell(Rc[i], scBoxMin_m,
                                                                 scBoxMax_m)) ? 1 : 0;
            charges[2 * part] += Qc[i];
            if (CoreFitSC::isSampled(IDc[i], coreStride[part])) {
                corePosition[part].push_back(CoreFitSC::toRestFrame(Rc[i], centre, gamma));
                coreCharge[part].push_back(Qc[i]);
                charges[2 * part + 1] += Qc[i];
            }
        }
        const bool sampled = (coreStride[0] > 1 || coreStride[1] > 1);
        if (sampled) {
            allreduce(charges, 4, std::plus<double>());
            ++ scStatistics_m.nearSampled;
        }
        std::vector<Vector_t> exact(flagged.size(), Vector_t(0.0));
        for (int part = 0; part < 2; ++part) {
            if (coreStride[part] > 1 && charges[2 * part + 1] != 0.0) {
                const double scale = charges[2 * part] / charges[2 * part + 1];
                for (double& q: coreCharge[part]) {
                    q *= scale;
                }
            }
            CoreFitSC::addDirectField(target, 0, flagged.size(), corePosition[part],
                                      coreCharge[part], 0.0, exact);
        }
        allreduce(&(exact[0](0)), 3 * static_cast<int>(flagged.size()), std::plus<double>());
        for (size_t k = 0; k < flagged.size(); ++k) {
            multipole[flagged[k]] = exact[k];
        }
        IpplTimings::stopTimer(exactCoreTimer_m);
    }

    // Lab-frame field of this rank's particles, with the coupling constant of the mesh
    // field; the MESH particles have their mesh field already.
    const double couplingConstant = getCouplingConstant();
    for (int j = 0; j < localCount; ++j) {
        const size_t row = offsets[myNode] + j;
        const size_t i = farIndex[j];
        if (particleClass[i] == Class_t::FAR) {
            Ef[i] = couplingConstant * CoreFitSC::toLabFrame(multipole[row] + farFar[row], gamma);
        } else {
            Ef[i] += couplingConstant * CoreFitSC::toLabFrame(farFar[row], gamma);
        }
    }

    if (farFarStride > 1) {
        ++ scStatistics_m.farFarSampled;
        std::ostringstream message;
        message << "the far-far field of the " << numLive << " far particles is summed over "
                << "about 1 in " << farFarStride << " of them, as " << numLive << " x "
                << numLive << " pairs exceed " << CoreFitSC::maxFarFarPairs;
        warnMeshFitOnce(FarFarWarning, message.str());
    }

    IpplTimings::stopTimer(farFieldTimer_m);
}


void PartBunch::checkUnmaskedInterpolation(const std::string& where) const {
    if (meshFit_m.type == MeshFitType::CORE && scStatistics_m.numOutside > 0) {
        throw GeneralClassicException(where,
                                      "FIELDSOLVER, MESHFIT=\"CORE\": " +
                                      std::to_string(scStatistics_m.numOutside) +
                                      " particles are outside the space-charge mesh, which "
                                      "the unmasked charge scatter and field gather cannot "
                                      "handle.");
    }
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
    INFOMSG(level3 << "BC set all periodic" << endl);
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
