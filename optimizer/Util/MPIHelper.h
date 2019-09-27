#ifndef __MPIHELPER_H__
#define __MPIHELPER_H__

#include <vector>
#include <map>
#include <sstream>
#include <mpi.h>

#include "Util/Types.h"


/// notify pilot about worker status
#define MPI_WORKER_STATUSUPDATE_TAG    0x11
/// notify pilot that work has been finished and results are ready to collect
#define MPI_WORKER_FINISHED_TAG        0x12
/// pilot notifies worker that he is ready to collect the results
#define MPI_WORKER_FINISHED_ACK_TAG    0x13

/// notifies coworkers about new jobs
#define MPI_COWORKER_NEW_JOB_TAG       0x14
#define MPI_WORKER_DIRECTORY_TAG       0x15
#define MPI_OPT_DIRECTORY_TAG          0x16


/// optimizer sends new job to pilot
#define MPI_OPT_NEW_JOB_TAG            0x21
/// pilot tells optimizer that results are ready to collect
#define MPI_OPT_JOB_FINISHED_TAG       0x22
/// optimizer notifies pilot that optimization has converged (EXIT)
#define MPI_OPT_CONVERGED_TAG          0x23


/// unique id of the job
#define MPI_WORK_JOBID_TAG             0x41
/// tags messages containing the size of results (requested vars)
#define MPI_WORK_SIZE_TAG              0x42
/// tags messages containing the size of the parameter list (for simulation)
#define MPI_WORK_SIZE_PARAMS           0x43

#define MPI_EXCHANGE_SOL_STATE_TAG          0x51
#define MPI_EXCHANGE_SOL_STATE_DATA_TAG     0x52
#define MPI_EXCHANGE_SOL_STATE_RES_SIZE_TAG 0x53
#define MPI_EXCHANGE_SOL_STATE_RES_TAG      0x54

/// global stop tag to exit poll loop (@see Poller)
#define MPI_STOP_TAG                   0x91

/// tag for exchanging serialized data
#define MPI_EXCHANGE_SERIALIZED_DATA_TAG        0x99


//FIXME: double information (use static vars or enums)
enum MPITag_t {
      WORKER_FINISHED_TAG = MPI_WORKER_FINISHED_TAG
    , OPT_NEW_JOB_TAG = MPI_OPT_NEW_JOB_TAG
    , OPT_CONVERGED_TAG = MPI_OPT_CONVERGED_TAG
    , WORKER_STATUSUPDATE_TAG = MPI_WORKER_STATUSUPDATE_TAG
    , REQUEST_FINISHED = MPI_OPT_JOB_FINISHED_TAG
    , EXCHANGE_SOL_STATE_TAG = MPI_EXCHANGE_SOL_STATE_TAG
    , EXCHANGE_SOL_STATE_RES_SIZE_TAG = MPI_EXCHANGE_SOL_STATE_RES_SIZE_TAG
};

/// Worker state is either idle or running
enum State_t {IDLE = 0, RUNNING = 1};

/// serializes params using Boost archive
void serialize(Param_t           params, std::ostringstream &os);
void serialize(reqVarContainer_t params, std::ostringstream &os);

/// deserializes params using Boost archive
void deserialize(char *buffer, Param_t           &params);
void deserialize(char *buffer, reqVarContainer_t &params);

/// broadcast params to all entities in comm
void MPI_Bcast_params(Param_t &params, size_t root, MPI_Comm comm);

/// broadcast requested variables to all entities in comm
void MPI_Bcast_reqvars(reqVarContainer_t reqvars, size_t root, MPI_Comm comm);


//FIXME: test
template<class Data_t>
void MPI_Send_serialized(Data_t data, size_t pid, MPI_Comm comm) {

    std::ostringstream os;
    serialize(data, os);
    size_t buf_size = os.str().length();

    MPI_Send(&buf_size, 1, MPI_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Send(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    delete[] buffer;
}

template<class Data_t>
void MPI_Recv_serialized(Data_t &data, size_t pid, MPI_Comm comm) {

    MPI_Status status;
    size_t buf_size = 0;
    MPI_Recv(&buf_size, 1, MPI_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    char *buffer = new char[buf_size];
    MPI_Recv(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    deserialize(buffer, data);
    delete[] buffer;
}



/**
 *  Serializes a parameter list and sends (MPI) to another processor.
 *
 *  @param params parameter list to serialize
 *  @param pid processor ID to send data to
 *  @param comm MPI communicator group used for sending the data
 */
void MPI_Send_params(Param_t params, size_t pid, MPI_Comm comm);

/**
 *  Serializes a parameter list and asynchronously sends (MPI) to another
 *  processor.
 *  The memory associated with returned buffers must be freed by the consumer.
 *
 *  @param params parameter list to serialize
 *  @param pid processor ID to send data to
 *  @param comm MPI communicator group used for sending the data
 *  @param req MPI request assigned with the ISend
 */
std::pair<size_t*, char*> MPI_ISend_params(Param_t params, size_t pid,
                                           MPI_Comm comm, MPI_Request *req);

/**
 *  Receives and unpacks a parameter list from another (MPI) processor.
 *
 *  @param params parameter list to store data
 *  @param pid processor ID to receive data from
 *  @param comm MPI communicator group used for receiving the data
 */
void MPI_Recv_params(Param_t &params, size_t pid, MPI_Comm comm);

/**
 *  Serializes requested variable list and sends (MPI) to another processor.
 *
 *  @param reqvars variable list list to serialize
 *  @param pid processor ID to send data to
 *  @param comm MPI communicator group used for sending the data
 */
void MPI_Send_reqvars(reqVarContainer_t reqvars, size_t pid, MPI_Comm comm);

/**
 *  Receives and unpacks a required variable list from another (MPI)
 *  processor.
 *
 *  @param reqvars variable list to store data
 *  @param pid processor ID to receive data from
 *  @param comm MPI communicator group used for receiving the data
 */
void MPI_Recv_reqvars(reqVarContainer_t &reqvars, size_t pid, MPI_Comm comm);

#endif
