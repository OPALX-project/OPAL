#ifndef __COMM_TOPOLOGY__
#define __COMM_TOPOLOGY__

#include <vector>
#include "mpi.h"


/// Specifies interface for topology policies
class CommTopology {

public:

    CommTopology(MPI_Comm comm = MPI_COMM_WORLD) : comm_(comm) {
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &num_procs_);
    }

    virtual ~CommTopology()
    {}

    /// every implementation must provide a discover method
    virtual void discover() = 0;

    int getRank() const { return rank_; }
    int getNP()   const { return num_procs_; }

    unsigned int getNumDimensions() const { return num_dims_; }
    unsigned int getCoreID()        const { return my_core_id_; }

    std::vector<unsigned int> getCoordinates() const { return coords_; }
    std::vector<unsigned int> getDimensions()  const { return dims_; }


private:

    MPI_Comm comm_;
    int rank_;
    int num_procs_;


protected:

    unsigned int hwID_;
    unsigned int num_dims_;
    unsigned int my_core_id_;

    std::vector<unsigned int> coords_;
    std::vector<unsigned int> dims_;

};

#endif
