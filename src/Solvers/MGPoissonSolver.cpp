
//#define DBG_STENCIL

#include "Solvers/MGPoissonSolver.h"

#include "Structure/BoundaryGeometry.h"
#include "ArbitraryDomain.h"
#include "EllipticDomain.h"
#include "BoxCornerDomain.h"

#include "Track/Track.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_Import.h"

#include "Teuchos_CommandLineProcessor.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

#include "Isorropia_Exception.hpp"
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraRedistributor.hpp"
#include "Isorropia_EpetraPartitioner.hpp"

#include <algorithm>

using Teuchos::RCP;
using Teuchos::rcp;
using Physics::c;

MGPoissonSolver::MGPoissonSolver ( PartBunch *beam,
                                   Mesh_t *mesh,
                                   FieldLayout_t *fl,
                                   std::vector<BoundaryGeometry *> geometries,
                                   std::string itsolver,
                                   std::string interpl,
                                   double tol,
                                   int maxiters,
                                   std::string precmode):
    geometries_m(geometries),
    tol_m(tol),
    maxiters_m(maxiters),
    Comm(Ippl::getComm()),
    itsBunch_m(beam),
    mesh_m(mesh),
    layout_m(fl) {

    domain_m = layout_m->getDomain();

    for (int i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        orig_nr_m[i] = domain_m[i].length();
    }

    if(itsolver == "CG") itsolver_m = AZ_cg;
    else if(itsolver == "BICGSTAB") itsolver_m = AZ_bicgstab;
    else if(itsolver == "GMRES") itsolver_m = AZ_gmres;
    else throw OpalException("MGPoissonSolver", "No valid iterative solver selected!");

    precmode_m = STD_PREC;
    if(precmode == "STD") precmode_m = STD_PREC;
    else if(precmode == "HIERARCHY") precmode_m = REUSE_HIERARCHY;
    else if(precmode == "REUSE") precmode_m = REUSE_PREC;
    else if(precmode == "NO") precmode_m = NO;

    tolerableIterationsCount_m = 2;
    numIter_m = -1;
    forcePreconditionerRecomputation_m = false;

    hasParallelDecompositionChanged_m = true;
    repartFreq_m = 1000;
    useRCB_m = false;
    if(Ippl::Info->getOutputLevel() > 3)
        verbose_m = true;
    else
        verbose_m = false;

    // Find CURRENT geometry
    currentGeometry = geometries_m[0];
    if(currentGeometry->getFilename() == "") {
        if(currentGeometry->getTopology() == "ELLIPTIC"){
            bp = new EllipticDomain(currentGeometry, orig_nr_m, hr_m, interpl);
        } else if (currentGeometry->getTopology() == "BOXCORNER") {
            bp = new BoxCornerDomain(currentGeometry->getA(), currentGeometry->getB(), currentGeometry->getC(), currentGeometry->getLength(),currentGeometry->getL1(), currentGeometry->getL2(), orig_nr_m, hr_m, interpl);
            bp->compute(itsBunch_m->get_hr());
        } else {
            throw OpalException("MGPoissonSolver::MGPoissonSolver",
                                "Geometry not known");
        }
    } else {
        NDIndex<3> localId = layout_m->getLocalNDIndex();
        if (localId[0].length() != domain_m[0].length() ||
            localId[1].length() != domain_m[1].length()) {
            throw OpalException("ArbitraryDomain::compute",
                                "The class ArbitraryDomain only works with parallelization\n"
                                "in z-direction.\n"
                                "Please set PARFFTX=FALSE, PARFFTY=FALSE, PARFFTT=TRUE in \n"
                                "the definition of the field solver in the input file.\n");
        }
        bp = new ArbitraryDomain(currentGeometry, orig_nr_m, hr_m, interpl);
    }

    Map = 0;
    A = Teuchos::null;
    LHS = Teuchos::null;
    RHS = Teuchos::null;
    MLPrec = 0;
    prec_m = Teuchos::null;

    numBlocks_m = Options::numBlocks;
    recycleBlocks_m = Options::recycleBlocks;
    nLHS_m = Options::nLHS;
    SetupMLList();
    SetupBelosList();
    problem_ptr = rcp(new Belos::LinearProblem<ST, MV, OP>);
    // setup Belos solver
    if(numBlocks_m == 0 || recycleBlocks_m == 0)
        solver_ptr = rcp(new Belos::BlockCGSolMgr<double, MV, OP>());
    else
        solver_ptr = rcp(new Belos::RCGSolMgr<double, MV, OP>());
    solver_ptr->setParameters(rcp(&belosList, false));
    convStatusTest = rcp(new Belos::StatusTestGenResNorm<ST, MV, OP> (tol));
    convStatusTest->defineScaleForm(Belos::NormOfRHS, Belos::TwoNorm);

    //all timers used here
    FunctionTimer1_m = IpplTimings::getTimer("BGF-IndexCoordMap");
    FunctionTimer2_m = IpplTimings::getTimer("computeMap");
    FunctionTimer3_m = IpplTimings::getTimer("IPPL to RHS");
    FunctionTimer4_m = IpplTimings::getTimer("ComputeStencil");
    FunctionTimer5_m = IpplTimings::getTimer("ML");
    FunctionTimer6_m = IpplTimings::getTimer("Setup");
    FunctionTimer7_m = IpplTimings::getTimer("CG");
    FunctionTimer8_m = IpplTimings::getTimer("LHS to IPPL");
}

void MGPoissonSolver::deletePtr() {
    delete Map;
    Map = nullptr;
    delete MLPrec;
    MLPrec = nullptr;
    A      = Teuchos::null;
    LHS    = Teuchos::null;
    RHS    = Teuchos::null;
    prec_m = Teuchos::null;
}

MGPoissonSolver::~MGPoissonSolver() {
    deletePtr ();
    solver_ptr = Teuchos::null;
    problem_ptr = Teuchos::null;
}

void MGPoissonSolver::computePotential(Field_t& /*rho*/, Vector_t /*hr*/, double /*zshift*/) {
    throw OpalException("MGPoissonSolver", "not implemented yet");
}

void MGPoissonSolver::computeMap(NDIndex<3> localId) {
    if (itsBunch_m->getLocalTrackStep()%repartFreq_m == 0){
        deletePtr();
        if(useRCB_m)
            redistributeWithRCB(localId);
        else
            IPPLToMap3D(localId);

        extrapolateLHS();
    }
}

void MGPoissonSolver::extrapolateLHS() {
// Aitken-Neville
// Pi0 (x) := yi , i = 0 : n
// Pik (x) := (x − xi ) Pi+1,k−1(x) − (x − xi+k ) Pi,k−1(x) /(xi+k − xi )
// k = 1, . . . , n, i = 0, . . . , n − k.
    //we also have to redistribute LHS

    if (Teuchos::is_null(LHS)){
        LHS = rcp(new Epetra_Vector(*Map));
        LHS->PutScalar(1.0);
    } else {
        RCP<Epetra_Vector> tmplhs = rcp(new Epetra_Vector(*Map));
        Epetra_Import importer(*Map, LHS->Map());
        tmplhs->Import(*LHS, importer, Add);
        LHS = tmplhs;
    }

    //...and all previously saved LHS
    std::deque< Epetra_Vector >::iterator it = OldLHS.begin();
    if(OldLHS.size() > 0) {
        int n = OldLHS.size();
        for(int i = 0; i < n; ++i) {
            Epetra_Vector tmplhs = Epetra_Vector(*Map);
            Epetra_Import importer(*Map, it->Map());
            tmplhs.Import(*it, importer, Add);
            *it = tmplhs;
            ++it;
        }
    }

    // extrapolate last OldLHS.size LHS to get a new start vector
    it = OldLHS.begin();
    if(nLHS_m == 0 || OldLHS.size()==0)
        LHS->PutScalar(1.0);
    else if(OldLHS.size() == 1)
        *LHS = *it;
    else if(OldLHS.size() == 2)
        LHS->Update(2.0, *it++, -1.0, *it, 0.0);
    else if(OldLHS.size() > 2) {
        int n = OldLHS.size();
        P = rcp(new Epetra_MultiVector(*Map, nLHS_m, false));
        for(int i = 0; i < n; ++i) {
           *(*P)(i) = *it++;
        }
        for(int k = 1; k < n; ++k) {
           for(int i = 0; i < n - k; ++i) {
              (*P)(i)->Update(-(i + 1) / (float)k, *(*P)(i + 1), (i + k + 1) / (float)k);
           }
        }
        *LHS = *(*P)(0);
     } else
        throw OpalException("MGPoissonSolver", "Invalid number of old LHS: " + std::to_string(OldLHS.size()));
}


// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation
// XXX: use matrix stencil in computation directly (no Epetra, define operators
// on IPPL GRID)
void MGPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {

    Inform msg("OPAL ", INFORM_ALL_NODES);
    nr_m[0] = orig_nr_m[0];
    nr_m[1] = orig_nr_m[1];
    nr_m[2] = orig_nr_m[2];

    bp->setGlobalMeanR(itsBunch_m->getGlobalMeanR());
    bp->setGlobalToLocalQuaternion(itsBunch_m->getGlobalToLocalQuaternion());
    bp->setNr(nr_m);

    NDIndex<3> localId = layout_m->getLocalNDIndex();

    IpplTimings::startTimer(FunctionTimer1_m);
    // Compute boundary intersections (only do on the first step)
    if(!itsBunch_m->getLocalTrackStep())
        bp->compute(hr, localId);
    IpplTimings::stopTimer(FunctionTimer1_m);

    // Define the Map
    INFOMSG(level3 << "* Computing Map..." << endl);
    IpplTimings::startTimer(FunctionTimer2_m);
    computeMap(localId);
    IpplTimings::stopTimer(FunctionTimer2_m);
    INFOMSG(level3 << "* Done." << endl);

    // Allocate the RHS with the new Epetra Map
    if (Teuchos::is_null(RHS))
        RHS = rcp(new Epetra_Vector(*Map));
    RHS->PutScalar(0.0);

    // // get charge densities from IPPL field and store in Epetra vector (RHS)
    if (verbose_m) {
        Ippl::Comm->barrier();
        msg << "* Node:" << Ippl::myNode() << ", Filling RHS..." << endl;
        Ippl::Comm->barrier();
    }
    IpplTimings::startTimer(FunctionTimer3_m);
    int id = 0;
    float scaleFactor = itsBunch_m->getdT();


    if (verbose_m) {
        msg << "* Node:" << Ippl::myNode() << ", Rho for final element: " << rho[localId[0].last()][localId[1].last()][localId[2].last()].get() << endl;

        Ippl::Comm->barrier();
        msg << "* Node:" << Ippl::myNode() << ", Local nx*ny*nz = " <<  localId[2].last() *  localId[0].last() *  localId[1].last() << endl;
        msg << "* Node:" << Ippl::myNode() << ", Number of reserved local elements in RHS: " << RHS->MyLength() << endl;
        msg << "* Node:" << Ippl::myNode() << ", Number of reserved global elements in RHS: " << RHS->GlobalLength() << endl;
        Ippl::Comm->barrier();
    }
    for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
                if (bp->isInside(idx, idy, idz))
                        RHS->Values()[id++] = 4*M_PI*rho[idx][idy][idz].get()/scaleFactor;
            }
        }
    }

    IpplTimings::stopTimer(FunctionTimer3_m);
    if (verbose_m) {
        Ippl::Comm->barrier();
        msg << "* Node:" << Ippl::myNode() << ", Number of Local Inside Points " << id << endl;
        msg << "* Node:" << Ippl::myNode() << ", Done." << endl;
        Ippl::Comm->barrier();
    }
    // build discretization matrix
    INFOMSG(level3 << "* Building Discretization Matrix..." << endl);
    IpplTimings::startTimer(FunctionTimer4_m);
    if(Teuchos::is_null(A))
        A = rcp(new Epetra_CrsMatrix(Copy, *Map,  7, true));
    ComputeStencil(hr, RHS);
    IpplTimings::stopTimer(FunctionTimer4_m);
    INFOMSG(level3 << "* Done." << endl);

#ifdef DBG_STENCIL
    EpetraExt::RowMatrixToMatlabFile("DiscrStencil.dat", *A);
#endif

    INFOMSG(level3 << "* Computing Preconditioner..." << endl);
    IpplTimings::startTimer(FunctionTimer5_m);
    if(!MLPrec) {
        MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m);
    } else if (precmode_m == REUSE_HIERARCHY) {
        MLPrec->ReComputePreconditioner();
    } else if (precmode_m == REUSE_PREC){
    }
    IpplTimings::stopTimer(FunctionTimer5_m);
    INFOMSG(level3 << "* Done." << endl);

    // setup preconditioned iterative solver
    // use old LHS solution as initial guess
    INFOMSG(level3 << "* Final Setup of Solver..." << endl);
    IpplTimings::startTimer(FunctionTimer6_m);
    problem_ptr->setOperator(A);
    problem_ptr->setLHS(LHS);
    problem_ptr->setRHS(RHS);
    if(Teuchos::is_null(prec_m))
        prec_m = Teuchos::rcp ( new Belos::EpetraPrecOp ( rcp(MLPrec,false)));
    problem_ptr->setLeftPrec(prec_m);
    solver_ptr->setProblem( problem_ptr);
    if (!problem_ptr->isProblemSet()) {
        if (problem_ptr->setProblem() == false) {
            ERRORMSG("Belos::LinearProblem failed to set up correctly!" << endl);
        }
    }
    IpplTimings::stopTimer(FunctionTimer6_m);
    INFOMSG(level3 << "* Done." << endl);

    double time = MPI_Wtime();

    INFOMSG(level3 << "* Solving for Space Charge..." << endl);
    IpplTimings::startTimer(FunctionTimer7_m);
    solver_ptr->solve();
    IpplTimings::stopTimer(FunctionTimer7_m);
    INFOMSG(level3 << "* Done." << endl);

    std::ofstream timings;
    if (true || verbose_m) {
        time = MPI_Wtime() - time;
        double minTime, maxTime, avgTime;
        Comm.MinAll(&time, &minTime, 1);
        Comm.MaxAll(&time, &maxTime, 1);
        Comm.SumAll(&time, &avgTime, 1);
        avgTime /= Comm.NumProc();
        if(Comm.MyPID() == 0) {
            char filename[50];
            sprintf(filename, "timing_MX%d_MY%d_MZ%d_nProc%d_recB%d_numB%d_nLHS%d", orig_nr_m[0], orig_nr_m[1], orig_nr_m[2], Comm.NumProc(), recycleBlocks_m, numBlocks_m, nLHS_m);
            timings.open(filename, std::ios::app);
            timings << solver_ptr->getNumIters() << "\t"
                    //<< time <<" "<<
                    << minTime << "\t"
                    << maxTime << "\t"
                    << avgTime << "\t"
                    << numBlocks_m << "\t"
                    << recycleBlocks_m << "\t"
                    << nLHS_m << "\t"
                    //<< OldLHS.size() <<"\t"<<
                    << std::endl;

            timings.close();
        }

    }
    // Store new LHS in OldLHS
    if(nLHS_m > 1) OldLHS.push_front(*(LHS.get()));
    if(OldLHS.size() > nLHS_m) OldLHS.pop_back();

    //now transfer solution back to IPPL grid
    IpplTimings::startTimer(FunctionTimer8_m);
    id = 0;
    rho = 0.0;
    for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
                  NDIndex<3> l(Index(idx, idx), Index(idy, idy), Index(idz, idz));
                  if (bp->isInside(idx, idy, idz))
                     rho.localElement(l) = LHS->Values()[id++]*scaleFactor;
            }
        }
    }
    IpplTimings::stopTimer(FunctionTimer8_m);

    if(itsBunch_m->getLocalTrackStep()+1 == (long long)Track::block->localTimeSteps.front()) {
        A = Teuchos::null;
        LHS = Teuchos::null;
        RHS = Teuchos::null;
        prec_m = Teuchos::null;
    }
}


void MGPoissonSolver::redistributeWithRCB(NDIndex<3> localId) {

    int numMyGridPoints = 0;

    for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
                if (bp->isInside(idx, idy, idz))
                    numMyGridPoints++;
             }
        }
     }

    Epetra_BlockMap bmap(-1, numMyGridPoints, 1, 0, Comm);
    Teuchos::RCP<const Epetra_MultiVector> coords = Teuchos::rcp(new Epetra_MultiVector(bmap, 3, false));

    double *v;
    int stride, stride2;

    coords->ExtractView(&v, &stride);
    stride2 = 2 * stride;

    for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
                if (bp->isInside(idx, idy, idz)) {
                    v[0] = (double)idx;
                    v[stride] = (double)idy;
                    v[stride2] = (double)idz;
                    v++;
                }
            }
        }
    }

    Teuchos::ParameterList paramlist;
    paramlist.set("Partitioning Method", "RCB");
    Teuchos::ParameterList &sublist = paramlist.sublist("ZOLTAN");
    sublist.set("RCB_RECTILINEAR_BLOCKS", "1");
    sublist.set("DEBUG_LEVEL", "1");

    Teuchos::RCP<Isorropia::Epetra::Partitioner> part = Teuchos::rcp(new Isorropia::Epetra::Partitioner(coords, paramlist));
    Isorropia::Epetra::Redistributor rd(part);
    Teuchos::RCP<Epetra_MultiVector> newcoords = rd.redistribute(*coords);

    newcoords->ExtractView(&v, &stride);
    stride2 = 2 * stride;
    numMyGridPoints = 0;
    std::vector<int> MyGlobalElements;
    for(int i = 0; i < newcoords->MyLength(); i++) {
        MyGlobalElements.push_back(bp->getIdx(v[0], v[stride], v[stride2]));
        v++;
        numMyGridPoints++;
    }

    Map = new Epetra_Map(-1, numMyGridPoints, &MyGlobalElements[0], 0, Comm);
}

void MGPoissonSolver::IPPLToMap3D(NDIndex<3> localId) {

    int NumMyElements = 0;
    std::vector<int> MyGlobalElements;

    for (int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        for (int idy = localId[1].first(); idy <= localId[1].last(); idy++) {
            for (int idx = localId[0].first(); idx <= localId[0].last(); idx++) {
                if (bp->isInside(idx, idy, idz)) {
                    MyGlobalElements.push_back(bp->getIdx(idx, idy, idz));
                    NumMyElements++;
                }
            }
        }
    }
    Map = new Epetra_Map(-1, NumMyElements, &MyGlobalElements[0], 0, Comm);
}

void MGPoissonSolver::ComputeStencil(Vector_t /*hr*/, Teuchos::RCP<Epetra_Vector> RHS) {

    A->PutScalar(0.0);

    int NumMyElements = Map->NumMyElements();
    int *MyGlobalElements = Map->MyGlobalElements();

    std::vector<double> Values(6);
    std::vector<int> Indices(6);

    for(int i = 0 ; i < NumMyElements ; i++) {

        int NumEntries = 0;

        double WV, EV, SV, NV, FV, BV, CV, scaleFactor=1.0;
        int W, E, S, N, F, B;

        bp->getBoundaryStencil(MyGlobalElements[i], WV, EV, SV, NV, FV, BV, CV, scaleFactor);
        RHS->Values()[i] *= scaleFactor;

        bp->getNeighbours(MyGlobalElements[i], W, E, S, N, F, B);
        if(E != -1) {
            Indices[NumEntries] = E;
            Values[NumEntries++] = EV;
        }
        if(W != -1) {
            Indices[NumEntries] = W;
            Values[NumEntries++] = WV;
        }
        if(S != -1) {
            Indices[NumEntries] = S;
            Values[NumEntries++] = SV;
        }
        if(N != -1) {
            Indices[NumEntries] = N;
            Values[NumEntries++] = NV;
        }
        if(F != -1) {
            Indices[NumEntries] = F;
            Values[NumEntries++] = FV;
        }
        if(B != -1) {
            Indices[NumEntries] = B;
            Values[NumEntries++] = BV;
        }

        // if matrix has already been filled (FillComplete()) we can only
        // replace entries

        if(A->Filled()) {
            // off-diagonal entries
            A->ReplaceGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
            // diagonal entry
            A->ReplaceGlobalValues(MyGlobalElements[i], 1, &CV, MyGlobalElements + i);
        } else {
            // off-diagonal entries
            A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
            // diagonal entry
            A->InsertGlobalValues(MyGlobalElements[i], 1, &CV, MyGlobalElements + i);
        }
    }

    A->FillComplete();

    A->OptimizeStorage();
}

void MGPoissonSolver::printLoadBalanceStats() {

    //compute some load balance statistics
    size_t myNumPart = Map->NumMyElements();
    size_t NumPart = Map->NumGlobalElements() * 1.0 / Comm.NumProc();
    double imbalance = 1.0;
    if(myNumPart >= NumPart)
        imbalance += (myNumPart - NumPart) / NumPart;
    else
        imbalance += (NumPart - myNumPart) / NumPart;

    double max = 0.0, min = 0.0, avg = 0.0;
    int minn = 0, maxn = 0;
    MPI_Reduce(&imbalance, &min, 1, MPI_DOUBLE, MPI_MIN, 0, Ippl::getComm());
    MPI_Reduce(&imbalance, &max, 1, MPI_DOUBLE, MPI_MAX, 0, Ippl::getComm());
    MPI_Reduce(&imbalance, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, Ippl::getComm());
    MPI_Reduce(&myNumPart, &minn, 1, MPI_INT, MPI_MIN, 0, Ippl::getComm());
    MPI_Reduce(&myNumPart, &maxn, 1, MPI_INT, MPI_MAX, 0, Ippl::getComm());

    avg /= Comm.NumProc();
    *gmsg << "LBAL min = " << min << ", max = " << max << ", avg = " << avg << endl;
    *gmsg << "min nr gridpoints = " << minn << ", max nr gridpoints = " << maxn << endl;
}

Inform &MGPoissonSolver::print(Inform &os) const {
    os << "* *************** M G P o i s s o n S o l v e r ************************************ " << endl;
    os << "* h " << hr_m << '\n';
    os << "* ********************************************************************************** " << endl;
    return os;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
