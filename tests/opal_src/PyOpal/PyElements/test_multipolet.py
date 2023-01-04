import math
import unittest 

import pyopal.elements.multipolet

class TestMultipoleT(unittest.TestCase):
    """
    Simple test on multipolet class
    """
    def make_multipolet(self):
        """Make a default straight multipole"""
        multipole = pyopal.elements.multipolet.MultipoleT()
        multipole.t_p = [0.5] # dipole
        multipole.left_fringe = 0.02
        multipole.right_fringe = 0.03
        multipole.length = 0.9
        multipole.horizontal_aperture = 0.15
        multipole.vertical_aperture = 0.25
        multipole.maximum_f_order = 3
        multipole.entrance_angle = 0.0 # not tested
        multipole.maximum_x_order = 3 # not tested
        multipole.variable_radius = 0 # not tested
        multipole.rotation = 0 # not tested
        multipole.angle = 0.0 # not tested
        multipole.bounding_box_length = 2.0 # doesn't appear to work
        multipole.delete_on_transvere_exit = True
        return multipole

    def test_on_axis_field(self):
        """Check the on axis field is reasonable"""
        multipole = self.make_multipolet()
        for zi, bz in [(-0.45, 0.25), (-0.0, 0.5), (0.45, 0.25)]:
            point = (0.0, 0.0, zi, 0.0)
            field = multipole.get_field_value(*point)
            self.assertAlmostEqual(field[2], bz)
            self.assertFalse(field[1])

    def test_bb(self):
        """Check the bounding box parsed correctly

        NOTE bb length does not appear to be used"""
        multipole = self.make_multipolet()
        plist = [
            ((0.07, 0.0, 0.0, 0.0), False),
            ((0.08, 0.0, 0.0, 0.0), True),
            ((0.0, 0.12, 0.0, 0.0), False),
            ((0.0, 0.13, 0.0, 0.0), True),
        ]
        for point, oob in plist:
            field = multipole.get_field_value(*point)
            self.assertEqual(oob, field[0], msg="Failed on point "+str(point)+" expected oob "+str(oob))
        field = multipole.get_field_value(0, 0, 10000, 0)
        self.assertEqual(False, field[0], msg="This test should have failed")
        multipole.delete_on_transvere_exit = False # disable bound checking
        field = multipole.get_field_value(100000, 100000, 10000, 100000)
        self.assertEqual(False, field[0], msg="Failed on point "+str(point)+" expected oob "+str(oob))

    def test_expansion(self):
        """Check the expansion parameters are parsed"""
        multipole = self.make_multipolet()
        multipole.maximum_f_order = 0
        multipole.maximum_x_order = 0 # this is related to opening angle stuff
        point = (0.05, 0.05, 0.45, 0.0)
        field_00 = multipole.get_field_value(*point)
        multipole.maximum_f_order = 1
        field_10 = multipole.get_field_value(*point)
        # just want to check they are different
        self.assertNotEqual(field_00[3], field_10[3])


if __name__ == "__main__":
    unittest.main()