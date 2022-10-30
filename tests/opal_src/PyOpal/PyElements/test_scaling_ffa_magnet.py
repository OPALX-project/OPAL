"""
Test ScalingFFAMagnet python implementation
"""

import math
import unittest
import pyopal.elements.scaling_ffa_magnet

class ScalingFFAMagnetTest(unittest.TestCase):
    """Test ScalingFFAMagnet"""
    def setUp(self):
        self.r0 = 4.0
        """Set some default magnet"""
        self.magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        self.magnet.b0 = 4.0
        self.magnet.r0 = self.r0
        self.magnet.field_index = 2
        self.magnet.tan_delta = math.tan(math.pi/4)
        self.magnet.max_vertical_power = 3
        #self.magnet.end_field_model =  
        self.magnet.end_length = 0.1
        self.magnet.centre_length = 1.1
        self.magnet.height = 0.61
        self.magnet.radial_neg_extent = 1.2
        self.magnet.radial_pos_extent = 1.3
        self.magnet.magnet_start = 1.5
        self.magnet.magnet_end = 2.5
        self.magnet.azimuthal_extent = 5.1 # that is extent from the centre to the end

    def get_cartesian_position(self, r, z, s):
        """
        Convert from magnet coordinates to cartesian coordinates. Magnet 
        coordinates are a circle, centred at x=-r0, having constant radius of 
        curvature r0. So x,y = 0,0 in cartesian coordinates correspond to 
        r,s=r0,0
        """
        theta = s/self.r0 # rad
        width = r*math.cos(theta)-self.r0
        longitudinal = r*math.sin(theta)
        return width, z, longitudinal, 0.0

    def test_get_field_value(self):
        """Check that we can get the field value okay"""
        for s in [0.1*i for i in range(75)]:
            point_cart = self.get_cartesian_position(self.r0, 0.0, s)
            oob, bx, by, bz, ex, ey, ez = self.magnet.get_field_value(*point_cart)
            #print("Tried", format(s, "6.4g"), "gave", oob, format(by, "6.4g"))

    def test_bounding_box_azimuthal(self):
        """Check the bounding box is okay in azimuthal coordinates"""
        self.magnet.tan_delta = 0.0
        for point_cyl in [ # just inside 8 corner points
            (self.r0-1.199, 0.304, 0.001), (self.r0-1.199, 0.304, 7.099), 
            (self.r0+1.299, 0.304, 0.001), (self.r0+1.299, 0.304, 7.099), 
            (self.r0-1.199, -0.304, 0.001), (self.r0-1.199, -0.304, 7.099), 
            (self.r0+1.299, -0.304, 0.001), (self.r0+1.299, -0.304, 7.099), 
            ]:
            point_cart = self.get_cartesian_position(*point_cyl)
            oob, bx, by, bz, ex, ey, ez = self.magnet.get_field_value(*point_cart)
            self.assertFalse(oob) # out of bounds
        for point_cyl in [ # take two opposite corner points and move off in each of three directions
            (self.r0-1.201, 0.304, 0.001), (self.r0-1.199, 0.306, 0.001), (self.r0-1.199, 0.306, -0.001), 
            (self.r0+1.301, -0.304, 7.199), (self.r0+1.299, -0.306, 7.199), (self.r0, -0.304, 7.201),
            ]:
            point_cart = self.get_cartesian_position(*point_cyl)
            oob, bx, by, bz, ex, ey, ez = self.magnet.get_field_value(*point_cart)
            self.assertTrue(oob, msg=f"Failed for {point_cyl}") # out of bounds


if __name__ == "__main__":
    unittest.main()
