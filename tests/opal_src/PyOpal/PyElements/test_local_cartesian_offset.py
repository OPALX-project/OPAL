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

"""Test LocalCartesianOffset python implementation"""
import unittest
import pyopal.objects.minimal_runner
import pyopal.elements.multipolet
import pyopal.elements.local_cartesian_offset
import pyopal.objects.field

class LocalCartesianOffsetTest(unittest.TestCase):
    """Tests for LocalCartesianOffset"""
    def test_init(self):
        """Check that we can build an offset okay"""
        offset = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        offset.end_position_x = 1.0
        self.assertEqual(offset.end_position_x, 1.0)
        offset.end_position_y = 2.0
        self.assertEqual(offset.end_position_y, 2.0)
        offset.end_normal_x = 3.0
        self.assertEqual(offset.end_normal_x, 3.0)
        offset.end_normal_y = 4.0
        self.assertEqual(offset.end_normal_y, 4.0)

    def test_placement(self):
        placement = OffsetPlacements()
        placement.offset = placement.null_drift()
        placement.execute_fork()

class OffsetPlacements(pyopal.objects.minimal_runner.MinimalRunner):
    """do placements"""
    def __init__(self):
        super().__init__()
        self.offset = None
        self.nominal_end_pos = (0, 1., 0, 0)
        self.nominal_end_dir = (0, 1, 0, 0)

    def vector_add(self, v1, v2):
        if len(v1) != len(v2):
            raise ValueError(f"Vector lengths of {v1} not equal to {v2}")
        v3 = tuple([v1_i+v2[i] for i, v1_i in enumerate(v1)])
        return v3



    def postprocess(self):
        pos1 = self.nominal_end_pos
        field = pyopal.objects.field.get_field_value(*pos1)
        print(pos1, field)

        pos2 = self.vector_add(self.nominal_end_pos, self.nominal_end_dir)
        field = pyopal.objects.field.get_field_value(*pos2)
        print(pos2, field)


    def make_element_iterable(self):
        return [self.offset, self.make_multipolet(1), self.make_multipolet(2)]

    def make_multipolet(self, bz):
        """Make a default straight multipole"""
        multipole = pyopal.elements.multipolet.MultipoleT()
        multipole.t_p = [bz] # dipole
        multipole.left_fringe = 0.0001
        multipole.right_fringe = 0.0001
        multipole.length = 0.001
        multipole.horizontal_aperture = 0.001
        multipole.vertical_aperture = 0.001
        multipole.maximum_f_order = 1
        multipole.entrance_angle = 0.0
        multipole.maximum_x_order = 1
        multipole.variable_radius = 0
        multipole.rotation = 0
        multipole.angle = 0.0
        multipole.bounding_box_length = 0.5
        multipole.delete_on_transvere_exit = False
        return multipole


if __name__ == "__main__":
    unittest.main()
