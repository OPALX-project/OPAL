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

import pyopal.objects.option
import pyopal.elements.probe


class ProbeRunner(pyopal.objects.minimal_runner.MinimalRunner):
    """Class to set up a probe and run it"""
    def __init__(self):
        """Set up the test"""
        super().__init__()
        self.probe = None
        self.options = None
        self.r0 = 0.9
        self.yprobe = 0.1*1e3 # that is in [mm]

    def make_options(self):
        """Don't want to force user to have hdf2py - so dump to ASCII"""
        self.options = pyopal.objects.option.Option()
        self.options.ascii_dump = True
        self.options.enable_hdf5 = False
        self.options.execute()

    def make_probe(self):
        """Make a probe, just offset from x axis so should catch first step"""
        self.probe = pyopal.elements.probe.Probe()
        self.probe.x_start = 0.0
        self.probe.x_end = self.r0*2*1000
        self.probe.y_start = self.yprobe
        self.probe.y_end = self.yprobe
        self.probe.output_filename = "test_probe"
        self.probe.set_opal_name("my_probe")

    def make_element_iterable(self):
        self.make_options()
        self.make_probe()
        return [self.probe]

class TestProbe(unittest.TestCase):
    def test_one(self):
        """The actual test - track through a probe and check we get an output"""
        probe_runner = ProbeRunner()
        probe_runner.execute_fork() # forked process for memory safety
        probe_file = os.path.join(probe_runner.tmp_dir, "test_probe.loss")
        with open(probe_file) as a_file:
            lines = [line for line in a_file.readlines()]
        self.assertEqual(len(lines), 3) # header and one track object
        words = lines[1].split()
        self.assertAlmostEqual(float(words[0]), probe_runner.r0)
        self.assertAlmostEqual(float(words[1]), probe_runner.yprobe*1e-3) # output is in [m]
        self.assertAlmostEqual(float(words[2]), 0.0)

if __name__ == "__main__":
    unittest.main()