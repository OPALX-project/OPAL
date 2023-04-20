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
import pyopal.elements.enge

class EngeTest(unittest.TestCase):
    """Test Enge"""
    def setUp(self):
        self.enge = pyopal.elements.enge.Enge()
        self.enge.x0 = 1.2
        self.enge.enge_lambda = 0.01
        self.enge.coefficients = [0.0, 0.3, 0.4, 0.5]

    def test_setup(self):
        """Check that we can get the field value okay"""
        self.assertEqual(self.enge.x0, 1.2)
        self.assertEqual(self.enge.enge_lambda, 0.01)
        self.assertEqual(self.enge.coefficients, [0.0, 0.3, 0.4, 0.5])

    def test_function(self):
        """Check that we can get values out of enge"""
        value = self.enge.function(0.0, 0)
        self.assertAlmostEqual(value, 1.0, 6)
        value = self.enge.function(1.2, 0)
        self.assertAlmostEqual(value, 0.5, 6)
        value = self.enge.function(-1.2, 0)
        self.assertAlmostEqual(value, 0.5, 6)
        self.enge.x0 = 2.3
        value = self.enge.function(2.3, 0)
        self.assertAlmostEqual(value, 0.5, 6)

    def test_name(self):
        self.enge.set_opal_name("bob")
        self.assertEqual(self.enge.get_opal_name(), "bob")

if __name__ == "__main__":
    unittest.main()
