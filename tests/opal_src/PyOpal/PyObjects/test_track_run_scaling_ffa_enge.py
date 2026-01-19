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
Module to test track run with scaling ffa lattice

This is fast enough to be a unit test; could go in the integration test 
framework but the integration test framework can't do it (yet).
"""

import os
import math
import unittest
import pyopal.elements.scaling_ffa_magnet
import pyopal.elements.asymmetric_enge
import pyopal.elements.enge
import pyopal.objects.field
import pyopal.objects.minimal_runner

C_LIGHT = 299792458 # m/s

class EngeScalingFFARunner(pyopal.objects.minimal_runner.MinimalRunner):
    """
    Module to test simulation of a scaling ffa (and by extension the entire 
    PyOpal workflow)
    """
    def __init__(self):
        """Here we set some basic lattice parameters"""
        super().__init__()
        self.r0 = 20.0
        self.run_name = "ScalingFFA"
        self.verbose = 0
        self.max_steps = 10
        self.mag_start = [0.5, 2.0, 3.0, 4.0]
        self.mag_end = [5.0, 6.0, 7.0, 8.0]
        self.mag_ctr = [0.3, 0.4, 0.5, 0.6]

    def make_element_iterable(self):
        """
        The lattice has two FFA magnets, a (F)ocusing magnet and a (D)efocusing
        magnet. In this example, the f magnet uses the default end field model,
        which is a tanh model. The d magnet uses a more general enge field
        model.
        """
        index = 0
        enge_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        enge_magnet.set_opal_name("enge_magnet")
        enge_magnet.set_attributes(
            r0=self.r0, b0=1.0, field_index=1, max_vertical_power=1,
            radial_neg_extent=0.1, radial_pos_extent=0.1, azimuthal_extent=1.0,
            end_length=99, centre_length=99,
            magnet_start=self.mag_start[index], magnet_end=self.mag_end[index],
            end_field_model="enge_end"
        )
        enge_end = pyopal.elements.enge.Enge()
        enge_end.set_opal_name("enge_end")
        enge_end.x0 = self.mag_ctr[index]
        enge_end.enge_lambda = 1e-6
        enge_end.coefficients = [0.0, 1.0]
        enge_end.update()
        enge_magnet.update_end_field()
        index += 1

        def_tanh_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        def_tanh_magnet.set_opal_name("def_tanh_magnet")
        def_tanh_magnet.set_attributes(
            r0=self.r0, b0=1.0, field_index=1, max_vertical_power=1,
            radial_neg_extent=0.1, radial_pos_extent=0.1, azimuthal_extent=1.0,
            end_length=1e-6, centre_length=self.mag_ctr[index],
            magnet_start=self.mag_start[index], magnet_end=self.mag_end[index],
        )
        index += 1

        tanh_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        tanh_magnet.set_opal_name("tanh_magnet")
        tanh_magnet.set_attributes(
            r0=self.r0, b0=1.0, field_index=1, max_vertical_power=1,
            radial_neg_extent=0.1, radial_pos_extent=0.1, azimuthal_extent=1.0,
            end_length=1e-6, centre_length=0.7,
            magnet_start=self.mag_start[index], magnet_end=self.mag_end[index],
            end_field_model="tanh_end"
        )
        enge_end = pyopal.elements.enge.Enge()
        enge_end.set_opal_name("tanh_end")
        enge_end.x0 = self.mag_ctr[index]
        enge_end.enge_lambda = 1e-6
        enge_end.coefficients = [0.0, 1.0]
        enge_end.update()
        enge_magnet.update_end_field()
        index += 1

        ass_enge_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        ass_enge_magnet.set_opal_name("ass_enge_magnet")
        ass_enge_magnet.set_attributes(
            r0=self.r0, b0=1.0, field_index=1, max_vertical_power=1,
            radial_neg_extent=0.1, radial_pos_extent=0.1, azimuthal_extent=1.0,
            end_length=99, centre_length=99,
            magnet_start=self.mag_start[index], magnet_end=self.mag_end[index],
            end_field_model="ass_enge_end"
        )
        ass_enge_end = pyopal.elements.asymmetric_enge.AsymmetricEnge()
        ass_enge_end.set_opal_name("ass_enge_end")
        ass_enge_end.x0_start = self.mag_ctr[index]-0.1
        ass_enge_end.x0_end = 0.1
        ass_enge_end.lambda_start = 1e-6
        ass_enge_end.lambda_end = 1e-6
        ass_enge_end.coefficients_start = [0.0, 1.0]
        ass_enge_end.coefficients_end = [0.0, 1.0]
        ass_enge_end.update()
        ass_enge_magnet.update_end_field()
        index += 1

        print(f"Built magnets with r0 {ass_enge_magnet.r0}")
        return [enge_magnet, def_tanh_magnet, tanh_magnet, ass_enge_magnet]

    def print_field(self):
        # magnet start at 1.0, 1.3, end at 2.0
        # magnet start at 3.5, 3.7 end at 4.5
        n_points = 500
        s_max = 25.0
        dx = 1e-3
        print (f"       s      phi        x        y field ....")
        for i in range(n_points+1):
            s = i/n_points*s_max
            field, phi, x, y = self.get_dipole_field(s)
            print (f"{s:8.4g} {phi:8.4g} {x:8.4g} {y:8.4g} {field[3]}")

    def get_dipole_field(self, s):
        phi = s/self.r0
        x = self.r0*math.cos(phi)
        y = self.r0*math.sin(phi)
        field = pyopal.objects.field.get_field_value(x, y, 0.0, 0.0)
        return field, phi, x, y

    def check_field(self):
        s0, ramping = 0.0, True
        for i in range(4):
            s0 += self.mag_start[i]
            for delta in [(-1e-3, 0.0), (0.0, 0.5), (1e-3, 1.0)]:
                s = s0+delta[0]
                field = self.get_dipole_field(s)[0][3]
                print(f"At s={s} m bz is {field} T")
                if abs(field - delta[1]) > 1e-3:
                    raise RuntimeError(f"Failed field check on magnet {i} ds: {delta[0]}")
            s0 += self.mag_end[i]-self.mag_start[i]

class TestTrackRun(unittest.TestCase):
    """Test class for track_run"""
    def test_run_one(self):
        """
        Test that we can run okay without an exception.

        If running from the command line, will spit out the OPAL log to screen.
        """
        runner = EngeScalingFFARunner()
        runner.postprocess = runner.print_field
        #runner.postprocess = runner.check_field
        runner.execute_fork()
        print("Finished\n\n")

if __name__ == "__main__":
    unittest.main()
