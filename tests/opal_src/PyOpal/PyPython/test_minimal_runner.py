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

import os
import unittest
import pyopal.objects.minimal_runner

class MinimalRunnerOverload(pyopal.objects.minimal_runner.MinimalRunner):
    """
    Test class that overloads a few of the functions in MinimalRunner for use in
    testing
    """
    def __init__(self):
        """Initialise"""
        super().__init__()
        self.preproc = os.path.join(self.tmp_dir, "preproc")
        self.postproc = os.path.join(self.tmp_dir, "postproc")

    def make_element_iterable(self):
        """This is just a dummy lattice, so just check that we can add
        some dummy drifts"""
        return [self.null_drift(), self.null_drift()]

    def preprocess(self):
        """Write a dummy  string at preprocess time"""
        # check that the element iterable was added to the line
        if len(self.line) != 4: # RingDef, minimal runner null drift, two more null drifts
            raise RuntimeError("Expected line of length 3, was actually "+str(len(self.line)))

        with open(self.preproc, "a", encoding="utf-8") as fout:
            fout.write("a")

    def postprocess(self):
        """Write a dummy  string at postprocess time"""
        with open(self.postproc, "a", encoding="utf-8") as fout:
            fout.write("b")

class MinimalRunnerTest(unittest.TestCase):
    def test_minimal_runner(self):
        """Test that minimal runner calls overloaded functions okay"""
        runner = MinimalRunnerOverload()
        # run 5 times to check that execution is suitably encapsulated
        for i in range(5):
            runner.execute_fork()
        # check that the runner preproc and postproc were called
        with open(runner.preproc, encoding="utf-8") as fin:
            self.assertEqual(len(fin.read()), 5)
        with open(runner.postproc, encoding="utf-8") as fin:
            self.assertEqual(len(fin.read()), 5)

if __name__ == "__main__":
    unittest.main()
