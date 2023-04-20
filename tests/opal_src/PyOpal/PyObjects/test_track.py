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

"""Test the track module"""
import unittest
import pyopal.objects.line
import pyopal.objects.beam
import pyopal.objects.track
import pyopal.objects.encapsulated_test_case

class TestTrack(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Quick check that we can run track okay"""

    def setUp(self):
        """Set up some data"""
        self.track = pyopal.objects.track.Track()

    def encapsulated_test_init(self):
        """Check I didn't make any typos"""
        self.track.line = "test_line"
        self.track.beam = "test_beam"
        self.track.dt = 1.0
        self.track.dt_space_charge = 2.0
        self.track.dtau = 3.0
        self.track.t0 = 4.0
        self.track.max_steps = [5]
        self.track.steps_per_turn = 6.0
        self.track.z_start = 7.0
        self.track.z_stop = [8.0]
        self.track.time_integrator = "integrator"
        self.track.map_order = 9

    def encapsulated_test_execute(self):
        """Check we can register the beam"""
        beam = pyopal.objects.beam.Beam()
        beam.set_opal_name("test_track::beam")
        beam.register()

        line = pyopal.objects.line.Line()
        line.set_opal_name("test_track::line")
        line.register()

        self.track.beam = "test_track::beam"
        self.track.line = "test_track::line"
        self.track.execute()


if __name__ == "__main__":
    unittest.main()
