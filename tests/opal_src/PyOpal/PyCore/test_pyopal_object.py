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
I picked on local_cartesian_offset as a simple element to test some generic
functionality for PyOpalObject (could do a mock-up, but this is simple)
"""

import unittest
import pyopal.elements.local_cartesian_offset

class PyOpalObjectTest(unittest.TestCase):
    """Test ScalingFFAMagnet"""
    def setUp(self):
        self.offset = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        self.offset.end_position_x = 1.0
        self.offset.end_position_y = 2.0
        self.offset.end_normal_x = 3.0
        self.offset.end_normal_y = 4.0

    def test_set_attributes(self):
        with self.assertRaises(RuntimeError): # called without class instance
            pyopal.elements.local_cartesian_offset.LocalCartesianOffset.set_attributes()
        with self.assertRaises(RuntimeError): # non-keyword argument
            self.offset.set_attributes(1)
        self.assertEqual(self.offset.end_position_x, 1.0)
        self.assertEqual(self.offset.end_position_y, 2.0)
        self.offset.set_attributes()
        self.assertEqual(self.offset.end_position_x, 1.0)
        self.assertEqual(self.offset.end_position_y, 2.0)
        self.offset.set_attributes(end_position_x = 20.0, end_position_y = 10.0)
        self.assertEqual(self.offset.end_position_x, 20.0)
        self.assertEqual(self.offset.end_position_y, 10.0)
        with self.assertRaises(RuntimeError): # not existing keyword argument
            self.offset.set_attributes(does_not_exist = 1.0)
        self.offset.set_attributes(end_position_x = 2.0, end_position_y = 1.0)
        self.assertEqual(self.offset.end_position_x, 2.0)
        self.assertEqual(self.offset.end_position_y, 1.0)

    def test_get_nonexistent_variable(self):
        with self.assertRaises(AttributeError):
            my_variable = self.offset.this_variable_does_not_exist_2

    def test_set_bad_type(self):
        try:
            self.offset.end_position_x = "not a float"
            self.assertTrue(False, msg="Should have failed")
        except Exception:
            pass

if __name__ == "__main__":
    unittest.main()
