/*
 *  Copyright (c) 2012, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OPALTESTUTILITIES_SILENCETEST_H_

#include <sstream>
#include <iostream>
#include "gtest/gtest.h"

namespace OpalTestUtilities {
    /** Shutup test output
     *
     *  If more than one is called, will shutup output on any alloc if it is loud
     *  and will make loud on any dealloc if it is quiet.
     */
    class FailureTester;

    class SilenceTest {
    public:
        SilenceTest();
        ~SilenceTest();

        void setFailed();
    private:
        SilenceTest(const SilenceTest& test); // disable default copy ctor

        std::ostringstream _debugOutput;
        static std::streambuf *_defaultCout;
        static std::streambuf *_defaultCerr;
        bool _failed;
        FailureTester *_failureTest;
    };

    class FailureTester: public ::testing::EmptyTestEventListener {
    public:
        FailureTester(SilenceTest *st):
            EmptyTestEventListener(),
            _test(st) { }

        ~FailureTester() { }

    private:
        FailureTester();

        virtual void OnTestPartResult(const ::testing::TestPartResult &test_part_result) {
            if (test_part_result.failed())
                _test->setFailed();
        }

        SilenceTest *_test;
    };
}

#endif //OPALTESTUTILITIES_SILENCETEST_H_