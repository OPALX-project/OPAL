import random
import unittest
from PyOpal.polynomial_map import PolynomialMap
from PyOpal.polynomial_coefficient import PolynomialCoefficient

class TestPolynomialMap(unittest.TestCase):
    def setUp(self):
        self.point_dim = 2
        self.poly_matrix = [[float(i) for i in range(0, 3)]+[0.],
                            [float(i) for i in range(3, 6)]+[0.],
                            [float(i) for i in range(6, 9)]+[0.],]
        self.ref_map = PolynomialMap(2, self.poly_matrix)

    def get_data(self, n_points):
        dimension = len(self.poly_matrix)
        points = [[random.randint(0, 10000)/10000. for i in range(self.point_dim)] \
                                                    for j in range(n_points)]
        values = [self.ref_map.evaluate(p) for p in points]
        return points, values

    def str_map(self, a_map):
        my_str = ""
        for row in a_map:
            for x in row:
                my_str += format(x, "8.2f") + " "
            my_str += "\n"
        return my_str

    def str_points_values(self, points, values):
        my_str = ""
        print 
        for i, p in enumerate(points):
            my_str += "p: "
            v = values[i]
            for x in p:
                my_str += format(x, "8.4g")+" "
            my_str += "   v: "
            for y in v:
                my_str += format(y, "8.4g")+" "
            my_str += "\n"

    def test_exact_solve(self):
        print "exact solve"
        points, values = self.get_data(4)
        test_map = PolynomialMap.exact_solve(points, values, 1)
        print "test"
        print self.str_map(test_map.get_coefficients_as_matrix())
        print "ref"
        print self.str_map(self.ref_map.get_coefficients_as_matrix())
        raise RuntimeError("This is failing")


    def test_least_squares(self):
        print "least squares"
        points, values = self.get_data(20)
        test_map = PolynomialMap.least_squares(points, values, 1)
        coeff_list = [
            PolynomialCoefficient
        ]
        print "test"
        print self.str_map(test_map.get_coefficients_as_matrix())
        print "ref"
        print self.str_map(self.ref_map.get_coefficients_as_matrix())


    def test_least_squares_coefficients(self):
        print "least squares coefficients"
        points, values = self.get_data(20)
        coeff_list = [
            PolynomialCoefficient([0, 1], i, 0.) for i in range(3)
        ]
        test_map = PolynomialMap.least_squares(points, values, 1, coeff_list)
        print "ref"
        print self.str_map(self.ref_map.get_coefficients_as_matrix())
        print "test 1"
        print self.str_map(test_map.get_coefficients_as_matrix())
        coeff_list += [
            PolynomialCoefficient([], 2, 0.)
        ]
        test_map = PolynomialMap.least_squares(points, values, 1, coeff_list)
        print "test 2"
        print self.str_map(test_map.get_coefficients_as_matrix())


    def test_everything(self):
        raise RuntimeError("More tests")

if __name__ == "__main__":
    unittest.main()