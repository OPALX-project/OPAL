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

"""Aim was to check that I pass processes argument to ippl.

As far as I can tell it is working okay (in that the arguments are passed to
ippl) but I can't see any effect on the execution. I do note that for my local
install, if I do something like 'mpirun -np 3 python test_globals.py --processes 4'
then ippl seems to hang."""

import unittest
import sys
import subprocess
import tempfile


def toy_import():
    """Import a toy simulation and run it. See if opal crashes"""
    import pyopal.objects.minimal_runner
    try:
        run = pyopal.objects.minimal_runner.MinimalRunner()
        run.verbose = 0
        run.execute()
    except RuntimeError:
        pass

class TestGlobals(unittest.TestCase):
    """Check that PyOpal globals set up works okay"""
    def run_one(self,  processes):
        """Run a toy simulation as a subproecss"""
        my_path = __file__
        args = ["python", my_path, "--processes", str(processes)]
        with tempfile.TemporaryFile() as stdout:
            my_proc = subprocess.Popen(args, stdout=stdout, stderr=subprocess.STDOUT)
            my_proc.wait()
            returnvalue = my_proc.returncode
            self.assertEqual(returnvalue, 0)

    def test_argparse(self):
        """Test that argument passing works"""
        self.run_one(3)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        toy_import()
    else:
        unittest.main()
