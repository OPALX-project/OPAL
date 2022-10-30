"""
Test AsymmetricEnge python implementation
"""

import math
import unittest
import pyopal.elements.asymmetric_enge

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
        value = self.enge.function(2.2825, 0)
        self.assertAlmostEqual(value, 0.5, 6)
        value = self.enge.function(-4.5, 0)
        self.assertAlmostEqual(value, 0.5, 6)

if __name__ == "__main__":
    unittest.main()
