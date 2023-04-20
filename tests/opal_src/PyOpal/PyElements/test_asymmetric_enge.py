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
Test AsymmetricEnge python implementation
"""

import math
import unittest
import pyopal.elements.asymmetric_enge
import pyopal.elements.enge

class AsymmetricEngeTest(unittest.TestCase):
    """Test AsymmetricEnge"""
    def setUp(self):
        self.enge = pyopal.elements.asymmetric_enge.AsymmetricEnge()
        self.enge.x0_start = 2.3
        self.enge.lambda_start = 0.02
        self.enge.coefficients_start = [0.0, 0.36, 0.46, 0.56]
        self.enge.x0_end = 1.2
        self.enge.lambda_end = 0.01
        self.enge.coefficients_end = [0.0, 0.3, 0.4, 0.5]

    def test_setup(self):
        """Check that we can get the field value okay"""
        self.assertEqual(self.enge.x0_start, 2.3)
        self.assertEqual(self.enge.lambda_start, 0.02)
        self.assertEqual(self.enge.coefficients_start, [0.0, 0.36, 0.46, 0.56])
        self.assertEqual(self.enge.x0_end, 1.2)
        self.assertEqual(self.enge.lambda_end, 0.01)
        self.assertEqual(self.enge.coefficients_end, [0.0, 0.3, 0.4, 0.5])

    def test_function(self):
        """Check that we can get values out of enge"""
        value = self.enge.function(0.0, 0)
        self.assertAlmostEqual(value, 1.0, 6)
        value = self.enge.function(2.3, 0)
        self.assertAlmostEqual(value, 0.5, 6)
        value = self.enge.function(-1.2, 0)
        self.assertAlmostEqual(value, 0.5, 6)

    def test_name(self):
        """Check that we can set and get name"""
        self.enge.set_opal_name("bob")
        self.assertEqual(self.enge.get_opal_name(), "bob")


if __name__ == "__main__":
    unittest.main()
