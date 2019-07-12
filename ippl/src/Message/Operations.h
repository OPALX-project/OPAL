#ifndef IPPL_MPI_OPERATIONS_H
#define IPPL_MPI_OPERATIONS_H

#include <functional>
#include <mpi.h>

template <class Op> MPI_Op get_mpi_op(Op op)
{
    return get_mpi_op( op );
}

#define IPPL_MPI_OP(CppOp, MPIOp)                   \
template <>                                         \
inline MPI_Op                                       \
get_mpi_op< CppOp >(CppOp) { return MPIOp; }

/* with C++14 we should be able
 * to simply write
 *
 * IPPL_MPI_OP(std::plus<>, MPI_SUM);
 *
 */
IPPL_MPI_OP(std::plus<double>, MPI_SUM);
// IPPL_MPI_OP(std::plus<long double>, MPI_SUM);
IPPL_MPI_OP(std::plus<int>, MPI_SUM);
IPPL_MPI_OP(std::plus<unsigned int>, MPI_SUM);
IPPL_MPI_OP(std::plus<size_t>, MPI_SUM);
IPPL_MPI_OP(std::plus<float>, MPI_SUM);
IPPL_MPI_OP(std::plus<long int>, MPI_SUM);


IPPL_MPI_OP(std::less<double>, MPI_MIN);
// IPPL_MPI_OP(std::less<long double>, MPI_MIN);
IPPL_MPI_OP(std::less<int>, MPI_MIN);
IPPL_MPI_OP(std::less<unsigned int>, MPI_MIN);
IPPL_MPI_OP(std::less<size_t>, MPI_MIN);
IPPL_MPI_OP(std::less<float>, MPI_MIN);
IPPL_MPI_OP(std::less<long int>, MPI_MIN);

IPPL_MPI_OP(std::greater<double>, MPI_MAX);
// IPPL_MPI_OP(std::greater<long double>, MPI_MAX);
IPPL_MPI_OP(std::greater<int>, MPI_MAX);
IPPL_MPI_OP(std::greater<unsigned int>, MPI_MAX);
IPPL_MPI_OP(std::greater<size_t>, MPI_MAX);
IPPL_MPI_OP(std::greater<float>, MPI_MAX);
IPPL_MPI_OP(std::greater<long int>, MPI_MAX);

IPPL_MPI_OP(std::logical_or<bool>, MPI_LOR);
IPPL_MPI_OP(std::logical_and<bool>, MPI_LAND);

#endif