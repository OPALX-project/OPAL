#include "gtest/gtest.h"
#include <gsl/gsl_errno.h>

#include "mpi.h"

#include "Utilities/OpalException.h"
#include "Utility/IpplInfo.h" // ippl

Ippl *ippl;
Inform* gmsg;

class NewLineAdder: public ::testing::EmptyTestEventListener {
    virtual void OnTestPartResult(const ::testing::TestPartResult &test_part_result) {
        if (test_part_result.failed())
            printf("\n");
    }
};

namespace {
    void errorHandlerGSL(const char *reason,
                         const char *file,
                         int /*line*/,
                         int /*gsl_errno*/) { // commented to stop gcc warning
        throw OpalException(file, reason);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    gmsg = new Inform("UnitTests: ", std::cerr);
    if (!gmsg) {
        return 1;
    }
    ippl = new Ippl(argc, argv);
    gsl_set_error_handler(&errorHandlerGSL);

    ::testing::TestEventListeners &listeners =
          ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new NewLineAdder);

    int test_out = RUN_ALL_TESTS();
    MPI_Finalize();

    return test_out;
}