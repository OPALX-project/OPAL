# Copyright (c) 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
#
# This file is part of OPAL.
#
# OPAL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Public License
# along with OPAL.  If not, see <https://www.gnu.org/licenses/>.

"""Module to run the unit tests"""

import tempfile
import inspect
import os
import sys
import unittest

def test_path():
    this_file = inspect.getfile(inspect.currentframe())
    this_dir = os.path.dirname(this_file)
    return os.path.join(this_dir)

class TestEnv():
    """Set up/tear down testing environment"""
    def __enter__(self):
        self.here = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)
        self.stream = open("test_output", "w+", encoding="utf-8")
        return self

    def __exit__(self, arg1, arg2, arg3):
        os.chdir(self.here)
        print(f"test output stored in dir {self.tmpdir}")

def main():
    """A very simple test runner script"""
    suite = unittest.defaultTestLoader.discover(
                start_dir = test_path(),
                pattern = "test*"
    )
    with TestEnv() as test_env:
        runner = unittest.TextTestRunner(stream=test_env.stream)
        result = runner.run(suite)
        print("Test output", result)
        test_env.stream.seek(0)
        print(test_env.stream.read())

    n_errors = len(result.errors)
    n_failures = len(result.failures)

    print(f"Ran tests with {n_errors} errors and {n_failures} failures")
    if n_errors or n_failures:
        print("Tests failed (don't forget to make install...)")
        sys.exit(1)
    else:
        print("Tests passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
