//
// Class MGPoissonSolver
//   This class contains methods for solving Poisson's equation for the
//   space charge portion of the calculation.
//
//   A smoothed aggregation based AMG preconditioned iterative solver for space charge
//   \see FFTPoissonSolver
//   \warning This solver is in an EXPERIMENTAL STAGE. For reliable simulations use the FFTPoissonSolver
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich,
//               2013 - 2015, Tülin Kaman, Paul Scherrer Institut, Villigen PSI, Switzerland,
//                      2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
//
// In 2020, the code was ported to the second generation Trilinos packages,
// i.e., Epetra --> Tpetra, ML --> MueLu. See also issue #507.
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
#ifndef MG_POISSON_SOLVER_H_
#define MG_POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"
#include "IrregularDomain.h"
//////////////////////////////////////////////////////////////

#include "mpi.h"

#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Teuchos_DefaultMpiComm.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>

#include "Teuchos_ParameterList.hpp"
#include "Algorithms/PartBunch.h"

typedef UniformCartesian<3, double> Mesh_t;
typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t> Field_t;

enum {
    STD_PREC,
    REUSE_PREC,
    REUSE_HIERARCHY
};

class BoundaryGeometry;

class MGPoissonSolver : public PoissonSolver {

public:
    typedef Tpetra::Vector<>                            TpetraVector_t;
    typedef Tpetra::MultiVector<>                       TpetraMultiVector_t;
    typedef Tpetra::Map<>                               TpetraMap_t;
    typedef Tpetra::Vector<>::scalar_type               TpetraScalar_t;
    typedef Tpetra::Vector<>::global_ordinal_type       TpetraGlobalOrdinal_t;
    typedef Tpetra::Operator<>                          TpetraOperator_t;
    typedef MueLu::TpetraOperator<>                     MueLuTpetraOperator_t;
    typedef Tpetra::CrsMatrix<>                         TpetraCrsMatrix_t;
    typedef Teuchos::MpiComm<int>                       Comm_t;

    typedef Teuchos::ParameterList                      ParameterList_t;

    typedef Belos::SolverManager<TpetraScalar_t,
                                 TpetraMultiVector_t,
                                 TpetraOperator_t>      SolverManager_t;

    typedef Belos::LinearProblem<TpetraScalar_t,
                                 TpetraMultiVector_t,
                                 TpetraOperator_t>      LinearProblem_t;



    MGPoissonSolver(PartBunch *beam,Mesh_t *mesh,
                    FieldLayout_t *fl,
                    std::vector<BoundaryGeometry *> geometries,
                    std::string itsolver, std::string interpl,
                    double tol, int maxiters, std::string precmode);

    ~MGPoissonSolver();

    /// given a charge-density field rho and a set of mesh spacings hr, compute
    /// the scalar potential in 'open space'
    /// \param rho (inout) scalar field of the potential
    /// \param hr mesh spacings in each direction
    void computePotential(Field_t &rho, Vector_t hr);
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    /// set a geometry
    void setGeometry(std::vector<BoundaryGeometry *> geometries);

    double getXRangeMin(unsigned short /*level*/) { return bp_m->getXRangeMin(); }
    double getXRangeMax(unsigned short /*level*/) { return bp_m->getXRangeMax(); }
    double getYRangeMin(unsigned short /*level*/) { return bp_m->getYRangeMin(); }
    double getYRangeMax(unsigned short /*level*/) { return bp_m->getYRangeMax(); }
    double getZRangeMin(unsigned short /*level*/) { return bp_m->getZRangeMin(); }
    double getZRangeMax(unsigned short /*level*/) { return bp_m->getZRangeMax(); }
    void test(PartBunchBase<double, 3>* /*bunch*/) { }
    /// useful load balance information
    void printLoadBalanceStats();

    void extrapolateLHS();

    Inform &print(Inform &os) const;

private:

    bool isMatrixfilled_m;

    // true if CG and GMRES; false if BiCGStab
    bool useLeftPrec_m;

    //TODO: we need to update this and maybe change attached
    //solver!
    /// holding the currently active geometry
    BoundaryGeometry *currentGeometry;

    /// container for multiple geometries
    std::vector<BoundaryGeometry *> geometries_m;

    /// flag notifying us that the geometry (discretization) has changed
    bool hasGeometryChanged_m;
    int repartFreq_m;
    /// flag specifying if we are verbose
    bool verbose_m;

    /// tolerance for the iterative solver
    double tol_m;
    /// maximal number of iterations for the iterative solver
    int maxiters_m;
    /// preconditioner mode
    int precmode_m;
    /// maximum number of blocks in Krylov space
    int numBlocks_m;
    /// number of vectors in recycle space
    int recycleBlocks_m;

    /// structure that holds boundary points
    std::unique_ptr<IrregularDomain> bp_m;

    /// right hand side of our problem
    Teuchos::RCP<TpetraVector_t> RHS;
    /// left hand side of the linear system of equations we solve
    Teuchos::RCP<TpetraVector_t> LHS;
    /// matrix used in the linear system of equations
    Teuchos::RCP<TpetraCrsMatrix_t> A;

    /// Map holding the processor distribution of data
    Teuchos::RCP<TpetraMap_t> map_p;
    /// communicator used by Trilinos
    Teuchos::RCP<const Comm_t> comm_mp;

    /// last N LHS's for extrapolating the new LHS as starting vector
    uint nLHS_m;
    Teuchos::RCP<TpetraMultiVector_t> P_mp;
    std::deque< TpetraVector_t > OldLHS;

    Teuchos::RCP<LinearProblem_t> problem_mp;
    Teuchos::RCP<SolverManager_t>  solver_mp;

    /// MueLu preconditioner object
    Teuchos::RCP<MueLuTpetraOperator_t> prec_mp;

    /// parameter list for the MueLu solver
    Teuchos::ParameterList MueLuList_m;
    /// parameter list for the iterative solver (Belos)
    Teuchos::ParameterList belosList;

    /// PartBunch object
    PartBunch *itsBunch_m;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;

    // domains for the various fields
    NDIndex<3> domain_m;

    /// mesh spacings in each direction
    Vector_t hr_m;
    /// current number of mesh points in each direction
    Vektor<int, 3> nr_m;
    /// global number of mesh points in each direction
    Vektor<int, 3> orig_nr_m;

    // timers
    IpplTimings::TimerRef FunctionTimer1_m;
    IpplTimings::TimerRef FunctionTimer2_m;
    IpplTimings::TimerRef FunctionTimer3_m;
    IpplTimings::TimerRef FunctionTimer4_m;
    IpplTimings::TimerRef FunctionTimer5_m;
    IpplTimings::TimerRef FunctionTimer6_m;
    IpplTimings::TimerRef FunctionTimer7_m;
    IpplTimings::TimerRef FunctionTimer8_m;

    void deletePtr();

    /// recomputes the map
    void computeMap(NDIndex<3> localId);

    /// converts IPPL grid to a 3D map
    /// \param localId local IPPL grid node indices
    void IPPLToMap3D(NDIndex<3> localId);

    /** returns a discretized stencil that has Neumann BC in z direction and
     * Dirichlet BC on the surface of a specified geometry
     * \param hr gridspacings in each direction
     * \param RHS right hand side might be scaled
     */
    void ComputeStencil(Vector_t hr, Teuchos::RCP<TpetraVector_t> RHS);



protected:
    /// Setup the parameters for the Belos iterative solver.
    void setupBelosList();

    /// Setup the parameters for the SAAMG preconditioner.
    void setupMueLuList();
};


inline Inform &operator<<(Inform &os, const MGPoissonSolver &fs) {
    return fs.print(os);
}

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:

