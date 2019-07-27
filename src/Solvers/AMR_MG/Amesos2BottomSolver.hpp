template <class Level>
Amesos2BottomSolver<Level>::Amesos2BottomSolver(std::string solvertype)
    : BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                Teuchos::RCP<amr::multivector_t>,
                                Level>()
    , solvertype_m(solvertype)
    , solver_mp(Teuchos::null)
{ }


template <class Level>
void Amesos2BottomSolver<Level>::solve(const Teuchos::RCP<mv_t>& x,
                                       const Teuchos::RCP<mv_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    solver_mp->solve(x.get(), b.get());
}


template <class Level>
void Amesos2BottomSolver<Level>::setOperator(const Teuchos::RCP<matrix_t>& A,
                                             Level* level_p)
{
    try {
        solver_mp = Amesos2::create<matrix_t, mv_t>(solvertype_m, A);
    } catch(const std::invalid_argument& ex) {
        *gmsg << ex.what() << endl;
    }
    
    solver_mp->symbolicFactorization();
    solver_mp->numericFactorization();

    this->isInitialized_m = true;
}


template <class Level>
std::size_t Amesos2BottomSolver<Level>::getNumIters() {
    return 1;   // direct solvers do only one step
}
