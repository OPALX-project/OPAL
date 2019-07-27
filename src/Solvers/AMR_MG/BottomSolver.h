#ifndef BOTTOM_SOLVER_H
#define BOTTOM_SOLVER_H

#include "AmrMultiGridDefs.h"

/// Abstract base class for all base level solvers
template <class Matrix, class Vector, class Level>
class BottomSolver {
    
public:

    BottomSolver() : isInitialized_m(false) { };
    
    virtual ~BottomSolver() { };
    
    /*!
     * Solves
     * \f[
     *      Ax = b
     * \f]
     * @param x left-hand side
     * @param b right-hand side
     */
    virtual void solve(const Vector& x,
                       const Vector& b) = 0;
    
    /*!
     * Set the system matrix
     * @param A system matrix
     */
    virtual void setOperator(const Matrix& A,
                             Level* level_p = nullptr) = 0;
    
    
    /*!
     * @returns the number of required iterations
     */
    virtual std::size_t getNumIters() = 0;

    bool hasOperator() const;


protected:
    bool isInitialized_m;
};


template <class Matrix, class Vector, class Level>
bool BottomSolver<Matrix, Vector, Level>::hasOperator() const {
    return isInitialized_m;
}

#endif
