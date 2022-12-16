"""
Test ScalingFFAMagnet python implementation
"""

import math
import unittest
import pyopal.elements.enge
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
        for s, by_ref in [(1.5, 2.0), (2.6, 2.0), (2.05, 4.0)]:
            point_cart = self.get_cartesian_position(self.r0, 0.0, s)
            oob, bx, by, bz, ex, ey, ez = self.magnet.get_field_value(*point_cart)
            #print("Tried", format(s, "6.4g"), "gave", oob, format(by, "6.4g"))
            self.assertAlmostEqual(by, by_ref, 3)

    def test_bounding_box_azimuthal(self):
        """Check the bounding box is okay in azimuthal coordinates"""
        self.magnet.tan_delta = 0.0
        for point_cyl in [ # just inside 8 corner points
            (self.r0-1.199, 0.304, 0.001), (self.r0-1.199, 0.304, 7.149), 
            (self.r0+1.299, 0.304, 0.001), (self.r0+1.299, 0.304, 7.149), 
            (self.r0-1.199, -0.304, 0.001), (self.r0-1.199, -0.304, 7.149), 
            (self.r0+1.299, -0.304, 0.001), (self.r0+1.299, -0.304, 7.149), 
            ]:
            point_cart = self.get_cartesian_position(*point_cyl)
            oob, bx, by, bz, ex, ey, ez = self.magnet.get_field_value(*point_cart)
            self.assertFalse(oob, msg=f"Failed for {point_cyl}") # out of bounds
        for point_cyl in [ # take two opposite corner points and move off in each of three directions
            (self.r0-1.201, 0.304, 0.001), (self.r0-1.199, 0.306, 0.001), (self.r0-1.199, 0.306, -0.001), 
            (self.r0+1.301, -0.304, 7.149), (self.r0+1.299, -0.306, 7.149), (self.r0+1.299, -0.304, 7.151),
            ]:
            point_cart = self.get_cartesian_position(*point_cyl)
            oob, bx, by, bz, ex, ey, ez = self.magnet.get_field_value(*point_cart)
            self.assertTrue(oob, msg=f"Failed for {point_cyl}") # out of bounds


    def test_end_field_model(self):
        """
        Check that the end field model is updated correctly. 

        Note I tried to set up some lazy/automatic updating for the endfieldmodel
        but I found it tricky to make it work. So in the end user has to cajole
        the end field model around the place. It's a bit clunky.
        """
        end_field = pyopal.elements.enge.Enge()
        end_field.enge_lambda = 1
        end_field.x0 = 2.5
        end_field.coefficients = [0.0, 1.0]
        end_field.set_opal_name("enge1")
        end_field.update()
        self.magnet.end_field_model = "enge1"
        point = self.get_cartesian_position(self.r0, 0.0, 4.6)
        by = self.magnet.get_field_value(*point)[2]
        # check normal operation
        self.assertAlmostEqual(by, 2.0, 2)

        # check updates in end field model propagate
        end_field.x0 = 2.2
        end_field.update()
        self.magnet.update_end_field()
        point = self.get_cartesian_position(self.r0, 0.0, 4.2)
        by = self.magnet.get_field_value(*point)[2]
        self.assertAlmostEqual(by, 2.0, 2)

        self.magnet.end_field_model = "engeNone"
        self.magnet.update_end_field()
        self.assertRaises(RuntimeError, self.magnet.get_field_value, *point)


if __name__ == "__main__":
    unittest.main()
