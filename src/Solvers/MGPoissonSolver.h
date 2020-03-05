//
// This class contains methods for solving Poisson's equation for the
// space charge portion of the calculation.
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#ifndef MG_POISSON_SOLVER_H_
#define MG_POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"
#include "IrregularDomain.h"
//////////////////////////////////////////////////////////////
#include "ml_include.h"

#include "mpi.h"

#if !defined(HAVE_ML_EPETRA)
#error Epetra support missing in Trilinos installation!
#endif

#if !defined(HAVE_ML_TEUCHOS)
#error Teuchos support missing in Trilinos installation!
#endif

#if !defined(HAVE_ML_AZTECOO)
#error Aztecoo support missing in Trilinos installation!
#endif

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_ParameterList.hpp"
#include "Algorithms/PartBunch.h"

#include "BelosTypes.hpp"

namespace Belos {
    template <class ScalarType, class MV, class OP>
    class LinearProblem;

    template <class ScalarType, class MV, class OP>
    class OperatorTraits;

    template<class ScalarType, class MV>
    class MultiVecTraits;

    template<class ScalarType, class MV, class OP>
    class SolverManager;

    class EpetraPrecOp;

    template <class ScalarType, class MV, class OP>
    class StatusTestGenResNorm;
}

namespace ML_Epetra {
    class MultiLevelPreconditioner;

    int SetDefaults(std::string ProblemType, Teuchos::ParameterList & List,
                    int * options, double * params, const bool OverWrite);
}


typedef UniformCartesian<3, double> Mesh_t;
typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t> Field_t;

enum {
    NO,
    STD_PREC,
    REUSE_PREC,
    REUSE_HIERARCHY
};

/**
 * \class MGPoissonSolver
 * \brief A smoothed aggregation based AMG preconditioned iterative solver for space charge
 * \see FFTPoissonSolver
 * \warning This solver is in an EXPERIMENTAL STAGE. For reliable simulations use the FFTPoissonSolver
 *
 */
class BoundaryGeometry;

class MGPoissonSolver : public PoissonSolver {

public:
    MGPoissonSolver(PartBunch *beam,Mesh_t *mesh, FieldLayout_t *fl, std::vector<BoundaryGeometry *> geometries, std::string itsolver, std::string interpl, double tol, int maxiters, std::string precmode);
    ~MGPoissonSolver();

    /// given a charge-density field rho and a set of mesh spacings hr, compute
    /// the scalar potential in 'open space'
    /// \param rho (inout) scalar field of the potential
    /// \param hr mesh spacings in each direction
    void computePotential(Field_t &rho, Vector_t hr);
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    /// set a geometry
    void setGeometry(std::vector<BoundaryGeometry *> geometries);

    /// force Solver to recompute Epetra_Map
    void recomputeMap() { hasParallelDecompositionChanged_m = true; }

    double getXRangeMin(unsigned short /*level*/) { return bp->getXRangeMin(); }
    double getXRangeMax(unsigned short /*level*/) { return bp->getXRangeMax(); }
    double getYRangeMin(unsigned short /*level*/) { return bp->getYRangeMin(); }
    double getYRangeMax(unsigned short /*level*/) { return bp->getYRangeMax(); }
    double getZRangeMin(unsigned short /*level*/) { return bp->getZRangeMin(); }
    double getZRangeMax(unsigned short /*level*/) { return bp->getZRangeMax(); }
    void test(PartBunchBase<double, 3>* /*bunch*/) { }
    /// useful load balance information
    void printLoadBalanceStats();

    void extrapolateLHS();

    Inform &print(Inform &os) const;

private:

    //TODO: we need to update this and maybe change attached
    //solver!
    /// holding the currently active geometry
    BoundaryGeometry *currentGeometry;

    /// container for multiple geometries
    std::vector<BoundaryGeometry *> geometries_m;

    /// flag notifying us that the geometry (discretization) has changed
    bool hasGeometryChanged_m;
    /// flag is set when OPAL changed decomposition of mesh
    bool hasParallelDecompositionChanged_m;
    int repartFreq_m;
    /// flag specifying if problem is redistributed with RCB
    bool useRCB_m;
    /// flag specifying if we are verbose
    bool verbose_m;

    /// tolerance for the iterative solver
    double tol_m;
    /// maximal number of iterations for the iterative solver
    int maxiters_m;
    /// iterative solver we are applying: CG, BiCGStab or GMRES
    int itsolver_m;
    /// preconditioner mode
    int precmode_m;
    /// number of iterations in the solve of the previous time step
    int numIter_m;
    /// percentage the iteration count can increase before recomputing the preconditioner
    int tolerableIterationsCount_m;
    /// force the solver to recompute the preconditioner
    bool forcePreconditionerRecomputation_m;
    /// maximum number of blocks in Krylov space
    int numBlocks_m;
    /// number of vectors in recycle space
    int recycleBlocks_m;

    /// structure that holds boundary points
    IrregularDomain *bp;

    /// right hand side of our problem
    Teuchos::RCP<Epetra_Vector> RHS;
    /// left hand side of the linear system of equations we solve
    Teuchos::RCP<Epetra_Vector> LHS;
    /// matrix used in the linear system of equations
    Teuchos::RCP<Epetra_CrsMatrix> A;

    /// ML preconditioner object
    ML_Epetra::MultiLevelPreconditioner *MLPrec;
    /// Epetra_Map holding the processor distribution of data
    Epetra_Map *Map;
    /// communicator used by Trilinos
    Epetra_MpiComm Comm;

    /// last N LHS's for extrapolating the new LHS as starting vector
    uint nLHS_m;
    Teuchos::RCP<Epetra_MultiVector> P;
    std::deque< Epetra_Vector > OldLHS;

    /// Solver (Belos RCG)
    typedef double                          ST;
    typedef Epetra_Operator                 OP;
    typedef Epetra_MultiVector              MV;
    typedef Belos::OperatorTraits<ST, MV, OP> OPT;
    typedef Belos::MultiVecTraits<ST, MV>    MVT;

    //Belos::LinearProblem<double, MV, OP> problem;
    typedef Belos::LinearProblem<ST, MV, OP> problem;
    Teuchos::RCP< problem > problem_ptr;

    typedef Belos::SolverManager<ST, MV, OP> solver;
    Teuchos::RCP< solver > solver_ptr;

    Teuchos::RCP< Belos::EpetraPrecOp > prec_m;
    Teuchos::RCP< Belos::StatusTestGenResNorm< ST, MV, OP > > convStatusTest;

    /// parameter list for the ML solver
    Teuchos::ParameterList MLList_m;
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

    /// recomputes the Epetra_Map
    void computeMap(NDIndex<3> localId);

    /// redistributes Map with RCB
    /// \param localId local IPPL grid node indices
    void redistributeWithRCB(NDIndex<3> localId);

    /// converts IPPL grid to a 3D Epetra_Map
    /// \param localId local IPPL grid node indices
    void IPPLToMap3D(NDIndex<3> localId);

    /** returns a discretized stencil that has Neumann BC in z direction and
     * Dirichlet BC on the surface of a specified geometry
     * \param hr gridspacings in each direction
     * \param Epetra_CrsMatrix holding the stencil
     * \param RHS right hand side might be scaled
     */
    void ComputeStencil(Vector_t hr, Teuchos::RCP<Epetra_Vector> RHS);



protected:

    /// Setup the parameters for the Belos iterative solver.
    inline void SetupBelosList() {
        belosList.set("Maximum Iterations", maxiters_m);
        belosList.set("Convergence Tolerance", tol_m);

        if(numBlocks_m != 0 && recycleBlocks_m != 0){//only set if solver==RCGSolMgr
            belosList.set("Num Blocks", numBlocks_m);               // Maximum number of blocks in Krylov space
            belosList.set("Num Recycled Blocks", recycleBlocks_m); // Number of vectors in recycle space
        }
        if(verbose_m) {
            belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails);
            belosList.set("Output Frequency", 1);
        } else
            belosList.set("Verbosity", Belos::Errors + Belos::Warnings);
      }

    /// Setup the parameters for the SAAMG preconditioner.
    inline void SetupMLList() {
        ML_Epetra::SetDefaults("SA", MLList_m, 0, 0, true);

        MLList_m.set("max levels", 8);
        MLList_m.set("increasing or decreasing", "increasing");

        // we use a V-cycle
        MLList_m.set("prec type", "MGV");

        // uncoupled aggregation is used (every processor aggregates
        // only local data)
        MLList_m.set("aggregation: type", "Uncoupled");

        // smoother related parameters
        MLList_m.set("smoother: type","Chebyshev");
        MLList_m.set("smoother: sweeps", 3);
        MLList_m.set("smoother: pre or post", "both");

        // on the coarsest level we solve with  Tim Davis' implementation of
        // Gilbert-Peierl's left-looking sparse partial pivoting algorithm,
        // with Eisenstat & Liu's symmetric pruning. Gilbert's version appears
        // as \c [L,U,P]=lu(A) in MATLAB. It doesn't exploit dense matrix
        // kernels, but it is the only sparse LU factorization algorithm known to be
        // asymptotically optimal, in the sense that it takes time proportional to the
        // number of floating-point operations.
        MLList_m.set("coarse: type", "Amesos-KLU");

        //XXX: or use Chebyshev coarse level solver
        // SEE PAPER FOR EVALUATION KLU vs. Chebyshev
        //MLList.set("coarse: sweeps", 10);
        //MLList.set("coarse: type", "Chebyshev");

        // Controls the amount of printed information.
        // Ranges from 0 to 10 (0 is no output, and
        // 10 is incredibly detailed output). Default: 0
        if(verbose_m)
            MLList_m.set("ML output", 10);

        // heuristic for max coarse size depending on number of processors
        int coarsest_size = std::max(Comm.NumProc() * 10, 1024);
        MLList_m.set("coarse: max size", coarsest_size);

    }

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

