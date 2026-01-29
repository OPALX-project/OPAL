#ifndef __NO_COMM_TOPOLOGY__
#define __NO_COMM_TOPOLOGY__

#include "mpi.h"
#include "Util/OptPilotException.h"
#include "Comm/Topology/CommTopology.h"

/// Simple policy when no topology is available or needed
class NoCommTopology : public CommTopology {

public:

    NoCommTopology(MPI_Comm comm = MPI_COMM_WORLD)
        : CommTopology(comm)
    {}

    virtual ~NoCommTopology()
    {}

    void discover() {
        throw OptPilotException("NoCommTopology::discoverTopology()",
                                "No topology policy selected!");
    }

};

#endif
