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

"""
Test Enge python implementation
"""

import math
import unittest
import pyopal.elements.polynomial_time_dependence

class PolynomialTimeDependenceTest(unittest.TestCase):
    """Test PolynomialTimeDependence"""
    def setUp(self):
        """Make a generic PTD"""
        self.ptd = pyopal.elements.polynomial_time_dependence.PolynomialTimeDependence()
        self.ptd.p0 = 1.0
        self.ptd.p1 = 2.0
        self.ptd.p2 = 3.0
        self.ptd.p3 = 4.0

    def test_setup(self):
        """Check that parameters are set okay"""
        self.assertEqual(self.ptd.p0, 1.0)
        self.assertEqual(self.ptd.p1, 2.0)
        self.assertEqual(self.ptd.p2, 3.0)
        self.assertEqual(self.ptd.p3, 4.0)

    def test_function(self):
        """Check that we can get correct values out"""
        for i in range(10):
            x = i*0.1
            value = self.ptd.function(x)
            test = sum([(i+1.0)*x**i for i in range(4)])
            self.assertAlmostEqual(test, value)

    def test_update(self):
        """No test - just check we don't get anything evil here"""
        self.ptd.update()

class PolynomialTimeDependenceTest2(unittest.TestCase):
    """Test PolynomialTimeDependence"""
    def setUp(self):
        """Make a generic PTD"""
        self.ptd = pyopal.elements.polynomial_time_dependence.PolynomialTimeDependence()
        self.ptd.coefficients = [1.0, 2.0, 3.0, 4.0]

    def test_setup(self):
        """Check that parameters are set okay"""
        self.assertEqual(self.ptd.coefficients, [1.0, 2.0, 3.0, 4.0])

    def test_function(self):
        """Check that we can get correct values out"""
        for i in range(10):
            x = i*0.1
            value = self.ptd.function(x)
            test = sum([(i+1.0)*x**i for i in range(4)])
            self.assertAlmostEqual(test, value)

    def test_update(self):
        """No test - just check we don't get anything evil here"""
        self.ptd.update()


if __name__ == "__main__":
    unittest.main()
