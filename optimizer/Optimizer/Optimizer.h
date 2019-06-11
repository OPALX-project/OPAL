#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include <vector>

#include "Util/Types.h"
#include "Pilot/Poller.h"

/**
 *  \class Optimizer
 *  \brief An abstract class defining the interface for all optimizer
 *         components.
 */
class Optimizer : protected Poller {

public:

    Optimizer(MPI_Comm comm) : Poller(comm) {}
    virtual ~Optimizer() {}

    /// type of bounds for design variables
    typedef std::vector< std::pair<double, double> > bounds_t;

    /// entry point for optimizer
    virtual void initialize() = 0;

protected:

    // propagate poller hooks
    virtual void setupPoll() = 0;
    virtual void prePoll() = 0;
    virtual void postPoll() = 0;
    virtual void onStop() = 0;
    virtual bool onMessage(MPI_Status status, size_t length) = 0;

private:

};

#endif
