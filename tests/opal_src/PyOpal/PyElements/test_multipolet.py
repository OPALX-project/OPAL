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

import math
import unittest
import sys
import os

has_ffa_field_mapper = False
try:
    import matplotlib
    import matplotlib.pyplot
    import pyopal.objects.ffa_field_mapper
    has_ffa_field_mapper = True
except ImportError:
    pass

import pyopal.objects.field
import pyopal.objects.minimal_runner
import pyopal.elements.multipolet
import pyopal.objects.field

class TestMultipoleTRunner(pyopal.objects.minimal_runner.MinimalRunner):
    """Test runner - I wanted to check the placement was okay"""
    def __init__(self):
        super().__init__()
        self.verbose = 0
        self.r0 = 2.0
        self.multipole = None
        self.length = 0.5
        self.bb_length = 0.6

    def make_element_iterable(self):
        """Place a 1.5 T dipole"""
        return [self.make_multipolet(1.5)]

    def make_multipolet(self, bz):
        """Make a default straight multipole"""
        multipole = pyopal.elements.multipolet.MultipoleT()
        multipole.set_attributes(
            t_p = [bz, 0.0, 0.0, 0.0/(0.05**3)], # dipole
            left_fringe = 0.01,
            right_fringe = 0.01,
            length = self.length,
            horizontal_aperture = 0.1,
            vertical_aperture = 0.1,
            maximum_f_order = 6,
            entrance_angle = 0.0,
            maximum_x_order = 6,
            variable_radius = 0,
            rotation = 0,
            angle = 0.0,
            bounding_box_length = self.bb_length,
            entry_offset = 0.0,
        )
        self.multipole = multipole
        return multipole

    @classmethod
    def test_field(cls, y_test, bz_ref):
        """Check that reference field is correct"""
        pos = (2.0, y_test, 0.0, 0.0) # x,y,z,t
        field = pyopal.objects.field.get_field_value(*pos)
        if abs(field[3]-bz_ref) > 1e-12:
            raise RuntimeError(f"test_multipole failed with pos {pos} field {field} bz_ref {bz_ref}")

    def postprocess(self):
        """Test placement - note that placement is now from entry point"""
        if has_ffa_field_mapper:
            self.plot()
        self.test_field(-0.301+0.25, 0.0)
        self.test_field(-0.299+0.25, 8.317278708416809e-05)
        self.test_field(0.299+0.25, 8.317278708416809e-05)
        self.test_field(0.301+0.25, 0.0)
        self.test_field(-0.25+0.25, 0.75)
        self.test_field(0.25+0.25, 0.75)

    def plot(self):
        mapper = pyopal.objects.ffa_field_mapper.FFAFieldMapper()
        mapper.x_points = [i*0.004+1.9 for i in range(50)]
        mapper.y_points = [i*0.01 for i in range(80)]
        mapper.z_point = 0.02
        mapper.field_map_cartesian(variable = "bz")
        mapper.field_map_cartesian(variable = "div_b")
        mapper.field_map_cartesian(variable = "curl_b")

    def _plot(self):
        """Plot the field"""
        if 'matplotlib' not in sys.modules:
            return
        x0 = 1.5
        x1 = 2.5
        y0 = 0.0
        y1 = 1.0
        n_x = 1000
        n_y = 1000
        x_list = []
        y_list = []
        bz_list = []
        dx = (x1-x0)/n_x
        dy = (y1-y0)/n_y
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        for i in range(n_x-1):
            for j in range(n_y-1):
                x = x0+(i+0.5)*dx
                y = y0+(j+0.5)*dy
                field = pyopal.objects.field.get_field_value(x, y, 0, 0)
                x_list.append(x)
                y_list.append(y)
                bz_list.append(field[3])
        bin_x = [x0+(x1-x0)*i/n_x for i in range(n_x)]
        bin_y = [y0+(y1-y0)*j/n_y for j in range(n_y)]
        axes.hist2d(x_list, y_list, weights=bz_list, bins=[bin_x, bin_y])
        fname = "test_multipolet_1.png"
        figure.savefig(fname)
        print(f"Saved figure in {os.getcwd()}/{fname}")

class TestMultipoleT(unittest.TestCase):
    """
    Simple test on multipolet class
    """
    def make_multipolet(self):
        """Make a default straight multipole"""
        multipole = pyopal.elements.multipolet.MultipoleT()
        multipole.t_p = [0.5] # 0.5 T dipole
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
        multipole.bounding_box_length = 2.0
        return multipole

    def test_on_axis_field(self):
        """Check the on axis field is reasonable"""
        multipole = self.make_multipolet()
        for zi, bz in [(0.0, 0.25), (0.45, 0.5), (0.9, 0.25)]:
            point = (0.0, 0.0, zi, 0.0)
            field = multipole.get_field_value(*point)
            self.assertAlmostEqual(field[2], bz)
            self.assertFalse(field[1])

    def test_bb(self):
        """Check the bounding box parsed correctly

        NOTE bb length does not appear to be used"""
        multipole = self.make_multipolet()
        plist = [ # should be rotated through 180 degrees?
            ((-0.07, 0.0, 0.45, 0.0), False),
            ((-0.08, 0.0, 0.45, 0.0), True),
            ((0.07, 0.0, 0.45, 0.0), False),
            ((0.08, 0.0, 0.45, 0.0), True),
            ((0.0, -0.12, 0.45, 0.0), False),
            ((0.0, -0.13, 0.45, 0.0), True),
            ((0.0, 0.12, 0.45, 0.0), False),
            ((0.0, 0.13, 0.45, 0.0), True),
            ((0.0, 0.0, -0.54, 0.0), False),
            ((0.0, 0.0, -0.56, 0.0), True),
            ((0.0, 0.0, 1.44, 0.0), False),
            ((0.0, 0.0, 1.46, 0.0), True),
        ]
        for point, oob in plist:
            field = multipole.get_field_value(*point)
            self.assertEqual(oob, field[0], msg="Failed on point "+str(point)+" expected oob "+str(oob))
        multipole.delete_on_transverse_exit = False # disable bound checking
        field = multipole.get_field_value(100000, 100000, 10000, 100000)
        self.assertEqual(False, field[0], msg="Failed on point "+str(point)+" expected oob "+str(oob))

    def test_expansion(self):
        """Check the expansion parameters are parsed"""
        multipole = self.make_multipolet()
        multipole.maximum_f_order = 1 # 0 is a bad value to use
        multipole.maximum_x_order = 0 # this is related to opening angle stuff
        point = (0.05, 0.05, 0.45, 0.0)
        field_00 = multipole.get_field_value(*point)
        multipole.maximum_f_order = 2
        field_10 = multipole.get_field_value(*point)
        # just want to check they are different
        self.assertNotEqual(field_00[3], field_10[3])

    def test_placement(self):
        """Check placement is okay"""
        runner = TestMultipoleTRunner()
        runner.execute_fork()

if __name__ == "__main__":
    unittest.main()