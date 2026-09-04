# Copyright (c) 2025, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
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

import math
import unittest
import pyopal.elements.sinusoidal_time_dependence

class SinusoidalTimeDependenceTest(unittest.TestCase):
    """Test SinusoidalTimeDependence"""
    def setUp(self):
        """Make a generic PTD"""
        self.ptd = pyopal.elements.sinusoidal_time_dependence.SinusoidalTimeDependence()
        self.ptd.frequencies = [8.0]
        self.ptd.amplitudes = [2.0]
        self.ptd.phase_offsets = [0.1]
        self.ptd.dc_offsets = [-1.0]

    def test_setup(self):
        """Check that parameters are set okay"""
        self.assertEqual(self.ptd.frequencies, [8.0])
        self.assertEqual(self.ptd.amplitudes, [2.0])
        self.assertEqual(self.ptd.phase_offsets, [0.1])
        self.assertEqual(self.ptd.dc_offsets, [-1.0])

    def test_function(self):
        """Check that we can get correct values out"""
        value = self.ptd.function(0.1)
        test = 2 / 2 * math.sin(2 * math.pi * 8 * 0.1 + 0.1) - 1
        self.assertAlmostEqual(test, value)

    def test_update(self):
        """No test - just check we don't get anything evil here"""
        self.ptd.update()


if __name__ == "__main__":
    unittest.main()
