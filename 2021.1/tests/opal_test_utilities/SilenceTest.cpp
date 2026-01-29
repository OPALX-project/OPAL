#include "opal_test_utilities/SilenceTest.h"
#include "Utility/IpplInfo.h"

std::streambuf *OpalTestUtilities::SilenceTest::_defaultCout = NULL;
std::streambuf *OpalTestUtilities::SilenceTest::_defaultCerr = NULL;

OpalTestUtilities::SilenceTest::SilenceTest():
    _failed(false) {
    IpplInfo::instantiateGlobals();
    if (_defaultCout == NULL ) {
        _defaultCout = std::cout.rdbuf();
        _defaultCerr = std::cerr.rdbuf();

        std::cout.rdbuf(_debugOutput.rdbuf());
        std::cerr.rdbuf(_debugOutput.rdbuf());

        ::testing::TestEventListeners& listeners =
              ::testing::UnitTest::GetInstance()->listeners();
        _failureTest = new FailureTester(this);
        listeners.Append(_failureTest);
    }
}

OpalTestUtilities::SilenceTest::~SilenceTest() { // return buffer to normal on delete
    if (_defaultCout != NULL) {
        std::cout.rdbuf(_defaultCout);
        std::cerr.rdbuf(_defaultCerr);
        _defaultCout = NULL;
        _defaultCerr = NULL;

        ::testing::TestEventListeners& listeners =
              ::testing::UnitTest::GetInstance()->listeners();
        listeners.Release(_failureTest);
        delete _failureTest;

        if (_failed) {
            std::cerr << _debugOutput.str() << std::endl;
        }
    }
}

void OpalTestUtilities::SilenceTest::setFailed() {
    _failed = true;
}