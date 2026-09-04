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

import unittest
import pyopal.elements.spline_time_dependence

class SplineTimeDependenceTest(unittest.TestCase):
    """Test SinusoidalTimeDependence"""
    def setUp(self):
        """Make a generic PTD"""
        self.ptd = pyopal.elements.spline_time_dependence.SplineTimeDependence()
        self.ptd.order = 1
        self.ptd.times = [0, 1]
        self.ptd.values = [0.1, 0.2]

    def test_setup(self):
        """Check that parameters are set okay"""
        self.assertEqual(self.ptd.order, 1)
        self.assertEqual(self.ptd.times, [0, 1])
        self.assertEqual(self.ptd.values, [0.1, 0.2])

    def test_function(self):
        """Check that we can get correct values out"""
        value = self.ptd.function(0)
        self.assertAlmostEqual(0.1, value)
        value = self.ptd.function(1)
        self.assertAlmostEqual(0.2, value)

    def test_update(self):
        """No test - just check we don't get anything evil here"""
        self.ptd.update()


if __name__ == "__main__":
    unittest.main()
