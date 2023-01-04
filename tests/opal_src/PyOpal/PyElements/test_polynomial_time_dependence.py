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


if __name__ == "__main__":
    unittest.main()
