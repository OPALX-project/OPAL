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

"""Test GlobalCartesianOffset python implementation"""
import unittest
import pyopal.objects.minimal_runner
import pyopal.elements.multipolet
import pyopal.elements.global_cartesian_offset
import pyopal.elements.local_cartesian_offset
import pyopal.objects.field

class GlobalCartesianOffsetTest(unittest.TestCase):
    """Tests for GlobalCartesianOffset"""
    def test_init(self):
        """Check that we can build an offset okay"""
        offset = pyopal.elements.global_cartesian_offset.GlobalCartesianOffset()
        offset.end_position_x = 1.0
        self.assertEqual(offset.end_position_x, 1.0)
        offset.end_position_y = 2.0
        self.assertEqual(offset.end_position_y, 2.0)
        offset.end_normal_x = 3.0
        self.assertEqual(offset.end_normal_x, 3.0)
        offset.end_normal_y = 4.0
        self.assertEqual(offset.end_normal_y, 4.0)

    def test_placement(self):
        """Check that placement works okay"""
        placement = OffsetPlacements()
        placement.offset = placement.null_drift()
        placement.execute_fork()

class OffsetPlacements(pyopal.objects.minimal_runner.MinimalRunner):
    """do placements"""
    def __init__(self):
        """initialisation"""
        super().__init__()
        self.verbose = 0

    def vector_add(self, vec1, vec2):
        """Add two vectors"""
        if len(vec1) != len(vec2):
            raise ValueError(f"Vector lengths of {vec1} not equal to {vec2}")
        vec3 = tuple([vec1_i+vec2[i] for i, vec1_i in enumerate(vec1)])
        return vec3

    def make_global_offset_1(self):
        """Make a global offset"""
        self.offset_1 = pyopal.elements.global_cartesian_offset.GlobalCartesianOffset()
        self.offset_1.set_opal_name("test_offset")
        self.offset_1.set_attributes(
            end_position_x = 1.0,
            end_position_y = 2.0,
            end_normal_x = 1/2**0.5,
            end_normal_y = 1/2**0.5
        )
        return self.offset_1

    def make_global_offset_2(self):
        """Make another global offset"""
        self.offset_2 = pyopal.elements.global_cartesian_offset.GlobalCartesianOffset()
        self.offset_2.set_opal_name("test_offset_2")
        self.offset_2.set_attributes(
            end_position_x = 0.5,
            end_position_y = 0.3,
            end_normal_x = 1.0,
            end_normal_y = 0.0
        )
        return self.offset_2


    def postprocess(self):
        """Check offsets are handled correctly"""
        tol = 1e-6
        # element 0 is a default "null" offset that comes with minimal_runner
        pos_offset_1 = pyopal.objects.field.get_element_end_position(1)
        assert abs(pos_offset_1[0] - self.offset_1.end_position_x) < tol
        assert abs(pos_offset_1[1] - self.offset_1.end_position_y) < tol
        norm_offset_1 = pyopal.objects.field.get_element_end_normal(1)
        assert abs(norm_offset_1[0] - self.offset_1.end_normal_x) < tol
        assert abs(norm_offset_1[1] - self.offset_1.end_normal_y) < tol
        # global offset 2 should place independently of offset 1
        pos_offset_2 = pyopal.objects.field.get_element_end_position(2)
        assert abs(pos_offset_2[0] - self.offset_2.end_position_x) < tol
        assert abs(pos_offset_2[1] - self.offset_2.end_position_y) < tol
        norm_offset_2 = pyopal.objects.field.get_element_end_normal(2)
        assert abs(norm_offset_2[0] - self.offset_2.end_normal_x) < tol
        assert abs(norm_offset_2[1] - self.offset_2.end_normal_y) < tol


    def make_element_iterable(self):
        """Set elements"""
        return [self.make_global_offset_1(), self.make_global_offset_2()]

if __name__ == "__main__":
    unittest.main()
