import unittest
from PyOpal.polynomial_coefficient import PolynomialCoefficient

class TestPolynomialCoefficient(unittest.TestCase):
    def setUp(self):
        pass
      
    def test_init(self):
        PolynomialCoefficient([0], 2, 3.2)

    def test_everything(self):
        raise RuntimeError("No tests")

if __name__ == "__main__":
    unittest.main()