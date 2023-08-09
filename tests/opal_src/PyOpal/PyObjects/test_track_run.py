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

"""Module to test track run (and by extension the entire PyOpal workflow)"""

import os
import unittest
import sys

import pyopal.elements.vertical_ffa_magnet
import pyopal.objects.minimal_runner
import pyopal.objects.encapsulated_test_case

class TrackRunner(pyopal.objects.minimal_runner.MinimalRunner):
    """
    Module to test track run, minimal track runner (and by extension the entire 
    PyOpal workflow)
    """
    def __init__(self):
        """Set up the test"""
        super().__init__()
        self.verbose = 0
        self.run_method = "CYCLOTRON-T"
        self.keep_alive = True
        self.beam_name = "DefaultBeam"
        self.distribution_name = ["DefaultDistribution"]
        self.field_solver_name = "DefaultFieldSolver"

    def make_track_run(self):
        run = pyopal.objects.track_run.TrackRun()
        run.method = self.run_method
        run.keep_alive = self.keep_alive
        run.beam_name = self.beam_name
        run.distribution = self.distribution_name
        run.field_solver = self.field_solver_name
        self.track_run = run

    def execute_fork(self):
        a_pid = os.fork()
        if a_pid == 0: # the child process
            try:
                self.execute()
            except RuntimeError:
                sys.excepthook(*sys.exc_info())
                os._exit(1)
            os._exit(0)
        else:
            retvalue = os.waitpid(a_pid, 0)[1]
        return retvalue



class TestTrackRun(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Test class for track_run"""
    def test_run_defaults(self):
        """
        Test that we can run okay without an exception.

        If running from the command line, will spit out the OPAL log to screen
        """
        track_runner = TrackRunner()
        self.assertEqual(track_runner.execute_fork(), 0) # Should throw??

    def test_catch_bad(self):
        track_runner = TrackRunner()
        self.assertEqual(track_runner.execute_fork(), 0)
        for var, substitute, retvalue in [
                ("run_method", "NO SUCH METHOD", 256),
                ("keep_alive", False, 0),
                ("beam_name", "CHEESE BEAM", 256),
                ("distribution_name", ["CHEESE DIST"], 256),
                ("field_solver_name", "BANANA", 256),
            ]:
            default = track_runner.__dict__[var]
            track_runner.__dict__[var] = substitute
            test_retvalue = track_runner.execute_fork()
            self.assertEqual(test_retvalue, retvalue)
            track_runner.__dict__[var] = default
        self.assertEqual(track_runner.execute_fork(), 0)

    def _encapsulated_test_run_two(self):
        """
        Test that we can run two track_runs without tripping on memory issues

        Test disabled - fails!
        """
        track_runner = TrackRunner()
        track_runner.execute()
        track_runner.make_track_run()
        track_runner.track_run.execute()

if __name__ == "__main__":
    unittest.main()
