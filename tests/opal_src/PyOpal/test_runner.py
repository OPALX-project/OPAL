"""Module to run the unit tests"""

import tempfile
import os
import sys
import unittest

class EncapsulatedRunner(unittest.TextTestRunner):
    """
    Copy of TextTestRunner but runs each test in it's own process
    """
    def __init__(self, stream=None):
        super().__init__(stream)

    def run(self, test):
        """
        Rather naive wrapper for TextTestRunner.run. Returns number of errors
        and failures in each child process.
        """
        print("Running", test)
        os.fork()
        a_pid = os.fork()
        if a_pid == 0: # the child process
            try:
                result = super().run(test)
                errors = len(result.errors)+len(result.failures)
                n_errors = len(result.errors)
                n_failures = len(result.failures)
            except:
                raise
            finally:
                os._exit(self.encode_return_value(n_errors, n_failures))
        else:
            print("starting", a_pid, test)
            return_value = os.waitpid(a_pid, 0)[1]
            n_errors, n_failures = self.decode_return_value(return_value)
            print("ending", a_pid, test)
        return n_errors, n_failures

    def encode_return_value(self, n_errors, n_failures):
        # retvalue is an 8 bit number; first 4 bits for n_errors, next 4 bits
        # for n_failures
        if n_errors > 15:
            n_errors = 15
        if n_failures > 15:
            n_failurse = 15
        return_value = n_errors + n_failures*16
        return return_value

    def decode_return_value(self, return_value):
        n_failures = return_value % 16
        n_errors = return_value - n_failures*16
        return n_failures, n_errors


def main():
    """A very simple test runner script"""
    suite = unittest.defaultTestLoader.discover(
                start_dir = "tests/opal_src/PyOpal/PyObjects/",
                pattern = "test*"
    )
    print("MAIN")
    with open("test", "w") as stream:
        runner = EncapsulatedRunner(stream=stream)
        n_errors, n_failures = runner.run(suite)
    with open("test") as stream:
        print(stream.read())
    print("MAIN 2")

    print(f"Ran tests with {n_errors} errors and {n_failures} failures")
    if n_errors or n_failures:
        print("Tests failed (don't forget to make install...)")
        sys.exit(n_errors+n_failures)
    else:
        print("Tests passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
