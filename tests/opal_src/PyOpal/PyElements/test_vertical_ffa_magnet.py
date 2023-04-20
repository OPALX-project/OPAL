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
Test VerticalFFAMagnet python implementation
"""

import math
import unittest
import pyopal.elements.vertical_ffa_magnet

class VerticalFFAMagnetTest(unittest.TestCase):
    """Test VerticalFFAMagnet"""
    def setUp(self):
        """Set some default magnet"""
        self.magnet = pyopal.elements.vertical_ffa_magnet.VerticalFFAMagnet()
        self.magnet.b0 = 4.0
        self.magnet.field_index = 2.0
        self.magnet.max_horizontal_power = 4
        self.magnet.centre_length = 1.0
        self.magnet.end_length = 0.01
        self.magnet.bb_length = 10.0
        self.magnet.height_neg_extent = 2.0
        self.magnet.height_pos_extent = 6.0
        self.magnet.width = 1.0

    def test_get_field_value(self):
        """Check that we can get the field value okay"""
        value = self.magnet.get_field_value(0.0, 0.0, 5.0, 0.0)
        self.assertEqual(len(value), 7)
        self.assertEqual(value[0], 0)
        self.assertAlmostEqual(value[1], 0)
        self.assertAlmostEqual(value[2], 4.0, 3)
        self.assertAlmostEqual(value[3], 0)
        self.assertAlmostEqual(value[4], 0)
        self.assertAlmostEqual(value[5], 0)
        self.assertAlmostEqual(value[6], 0)

        value = self.magnet.get_field_value(0.0, 1.0, 5.0, 0.0)
        self.assertAlmostEqual(value[2], 4.0*math.exp(2.0*1.0), 3)
        value = self.magnet.get_field_value(0.0, 0.0, 4.5, 0.0)
        self.assertAlmostEqual(value[2], 2.0, 3)
        value = self.magnet.get_field_value(0.0, 0.0, 5.5, 0.0)
        self.assertAlmostEqual(value[2], 2.0, 3)
        value = self.magnet.get_field_value(0.2, 0.0, 4.5, 0.0)
        self.assertAlmostEqual(value[2], 1.84, 2)

    def test_bounding_box(self):
        """Check that we can set up the bounding box okay"""
        for point in [(0.0, 0.0, -0.01, 0.0), (0.0, 0.0, 10.01, 0.0),
                      (0.0, -2.01, 5.0, 0.0), (0.0, 6.01, 5.0, 0.0),
                      (-0.51, 5.0, 0.0, 0.0), (0.51, 5.0, 0.0, 0.0),]:
            value = self.magnet.get_field_value(*point)
            self.assertTrue(value[0], msg=f"{point} {value}")
        for point in [(0.0, 0.0, 0.01, 0.0), (0.0, 0.0, 9.99, 0.0),
                      (0.0, -1.99, 5.0, 0.0), (0.0, 5.99, 5.0, 0.0),
                      (-0.49, 5.0, 0.0, 0.0), (0.49, 5.0, 0.0, 0.0),]:
            value = self.magnet.get_field_value(*point)
            self.assertFalse(value[0], msg=f"{point} {value}")

if __name__ == "__main__":
    unittest.main()
