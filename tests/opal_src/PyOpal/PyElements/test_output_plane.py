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

import os
import unittest
import math
import pyopal.objects.minimal_runner

import pyopal.objects.field
import pyopal.objects.option
import pyopal.elements.output_plane
import pyopal.elements.multipolet
import pyopal.elements.global_cartesian_offset
import pyopal.elements.polynomial_time_dependence
import pyopal.elements.variable_rf_cavity


class OutputPlaneRunner(pyopal.objects.minimal_runner.MinimalRunner):
    """Class to set up a probe and run it"""
    def __init__(self, number_of_steps=100, beta_rel=0.001, bz=0.0, algorithm="RK4", verbose=0):
        """Set up the test"""
        super().__init__()
        self.verbose = verbose
        self.probe = None
        self.options = None
        self.r0 = 0.9
        self.yprobe = 0.1 # that should be in [m]
        self.algorithm = algorithm

        self.bz = bz
        # Tracking is set up so the beam will go 0.11 metres
        self.steps_per_turn = number_of_steps-1
        self.max_steps = number_of_steps
        self.time_per_turn = 0.11/(beta_rel*3e8) # [seconds] beam goes 0.11 metres
        self.momentum = beta_rel/(1-beta_rel**2)*self.mass # [GeV/c]
        self.distribution_str = """1
0.0 0.0 0.0 0.0 0.0 0.0
"""
        self.setup_plane()

    def make_options(self):
        """Don't want to force user to have hdf2py - so dump to ASCII"""
        self.options = pyopal.objects.option.Option()
        self.options.ascii_dump = True
        self.options.enable_hdf5 = False
        self.options.execute()

    def setup_plane(self):
        """Make a probe, just offset from x axis so should catch first step"""
        self.output_plane = pyopal.elements.output_plane.OutputPlane()
        self.output_plane.set_attributes(
            x_start = 0.1,
            x_end = self.r0*2,
            y_start = self.yprobe,
            y_end = self.yprobe,
            placement_style = "PROBE",
            output_filename = "test_plane",
            algorithm = self.algorithm,
            verbose_level = self.verbose
        )
        self.output_plane.set_opal_name("my_plane")

    def make_dipole(self):
        dipole = pyopal.elements.multipolet.MultipoleT()
        dipole.set_attributes(
            t_p = [self.bz, self.bz*100],
            left_fringe = 0.01,
            right_fringe = 0.01,
            length = 0.5,
            horizontal_aperture = 100,
            vertical_aperture = 100,
            bounding_box_length = 1.0
        )
        return [dipole]

    @classmethod
    def make_time_dependence(cls, name, pol0, pol1):
        """Make a time dependence"""
        ptd = pyopal.elements.polynomial_time_dependence.PolynomialTimeDependence()
        ptd.p0 = pol0
        ptd.p1 = pol1
        ptd.set_opal_name(name)
        ptd.update()
        return ptd

    def make_rf(self):
        rf = pyopal.elements.variable_rf_cavity.VariableRFCavity()
        phase = self.make_time_dependence("phase", math.pi/2, 0.0)
        # nb: 4 GeV protons; significant energy change in 0.1 m -> O(GV/m)
        voltage = self.make_time_dependence("voltage", self.bz*1000.0, 0.0)
        # time step is ~ ns so frequency ~ GHz is appropriate
        frequency = self.make_time_dependence("frequency", 1e3, 0.0)
        rf.set_attributes(
            phase_model = "phase",
            amplitude_model = "voltage",
            frequency_model = "frequency",
            height = 10.0,
            width = 10.0,
            length = 100.0,
        )
        offset = pyopal.elements.global_cartesian_offset.GlobalCartesianOffset()
        offset.set_attributes(
            end_position_x = 0.9,
            end_position_y = -0.5,
            end_normal_x = 0.0,
            end_normal_y = 1.0,
        )
        return [offset, rf]

    def make_element_iterable(self):
        self.make_options()
        return [self.output_plane]+self.make_dipole()+self.make_rf()

class TestOutputPlane(unittest.TestCase):

    def run_one(self, n_steps, beta_rel, bz, algorithm):
        plane_runner = OutputPlaneRunner(n_steps, beta_rel, bz, algorithm, self.verbose)
        plane_runner.execute_fork() # forked process for memory safety
        return self.read_lines(plane_runner)

    def read_lines(self, plane_runner, name = "test_plane"):
        plane_file = os.path.join(plane_runner.tmp_dir, name+".loss")
        try:
            with open(plane_file) as a_file:
                lines = [line for line in a_file.readlines()]
        except FileNotFoundError:
            print(f"File {plane_file} not found - ls returns:")
            print("   ", os.listdir(plane_runner.tmp_dir))
            raise
        return lines

    def test_no_field_rk4(self):
        """Simple test - track through a probe and check we get an output"""
        lines = self.run_one(1000, 0.001, 0.0, "INTERPOLATION")
        self.assertEqual(len(lines), 2) # header and one track object
        words = lines[1].split()
        self.assertAlmostEqual(float(words[0]), 0.9)
        self.assertAlmostEqual(float(words[1]), 0.1) # output is in [m]
        self.assertAlmostEqual(float(words[2]), 0.0)
        self.assertAlmostEqual(float(words[3]), 0.0)
        self.assertAlmostEqual(float(words[4]), 0.001) # beta_gamma
        self.assertAlmostEqual(float(words[5]), 0.0)

    def test_field_rk4(self):
        """
        Track through a probe and check we get an output; check that RK4 is
        consistent with linear interpolation in the limit of many interpolation
        steps. Note fields are set up just to be "reasonably" complicated i.e. I
        want to check fields in several different directions, time varying, etc
        all works okay.
        """
        int_lines = self.run_one(1000, 0.9, 2.0, "INTERPOLATION")
        rk4_lines = self.run_one(3, 0.9, 2.0, "RK4")
        if self.verbose:
            print(rk4_lines)
            print("Interpolation\n   ", int_lines[1])
            print("RK4\n   ", rk4_lines[1])
        self.assertEqual(len(int_lines), 2) # header and one track object
        self.assertEqual(len(rk4_lines), 2) # header and one track object
        rk4_words = rk4_lines[1].split()
        int_words = int_lines[1].split()
        self.assertAlmostEqual(float(rk4_words[3]), float(int_words[3]), 4)
        self.assertAlmostEqual(float(rk4_words[4]), float(int_words[4]), 4) # beta_gamma
        self.assertAlmostEqual(float(rk4_words[5]), float(int_words[5]), 4)
        self.assertAlmostEqual(float(rk4_words[0]), float(int_words[0]), 9)
        self.assertAlmostEqual(float(rk4_words[1]), float(int_words[1]), 9) # output is in [m]
        self.assertAlmostEqual(float(rk4_words[2]), float(int_words[2]), 9)
        self.assertAlmostEqual(float(rk4_words[9])*1e9, float(int_words[9])*1e9, 3) # ns

    def test_reference_alignment_particle(self):
        """Check plane realignment works"""
        plane_runner = OutputPlaneRunner(100, 0.001, 0.0, "INTERPOLATION", self.verbose)
        plane_runner.output_plane.reference_alignment_particle = 1
        plane_runner.output_plane.tolerance = 1e-10
        # momentum units are eV(!)
        dpx = 0.0001*0.938272e9 # dxds = 0.1
        plane_runner.distribution_str = """3
0.0 0.0 0.0 0.0 0.0 0.0
0.0 """+str(dpx)+""" 0.0 0.0 0.0 0.0
-0.01 0.0 0.0 0.0 0.0 0.0
"""
        plane_runner.execute_fork() # forked process for memory safety
        lines = self.read_lines(plane_runner)
        if self.verbose:
            print("Reference plane test")
            for line in lines:
                print(line[:-1])
        self.assertEqual(len(lines), 4)
        dxds = 0.1 # the alignment particle angle
        dx = 0.02 # alignment particle offset
        test_particle = lines[3].split()
        self.assertAlmostEqual(float(test_particle[0]), 0.89)
        self.assertAlmostEqual(float(test_particle[1])-0.1, dx*dxds)

    def test_rectangular_aperture(self):
        """Check we cut approporiately with rectangular aperture"""
        width = 0.01
        plane_runner = OutputPlaneRunner(100, 0.001, 0.0, "INTERPOLATION", self.verbose)
        plane_runner.output_plane.placement_style = "CENTRE-NORMAL"
        plane_runner.output_plane.centre = [0.9, 0.1, 0.0]
        plane_runner.output_plane.normal = [1.0, 1.0, 0.0]
        plane_runner.output_plane.height = 0.1
        plane_runner.output_plane.width = width*2
        plane_runner.distribution_str = f"""9
0.0 0.0 0.0 0.0 0.0 0.0
{(width-0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
{(width+0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
{(-width-0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
{(-width+0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.049 0.0
0.0 0.0 0.0 0.0 0.051 0.0
0.0 0.0 0.0 0.0 -0.049 0.0
0.0 0.0 0.0 0.0 -0.051 0.0
"""

        plane_runner.execute_fork() # forked process for memory safety
        lines = self.read_lines(plane_runner)
        if self.verbose:
            print("Reference plane test")
            for line in lines:
                print(line[:-1])
        self.assertEqual(len(lines), 6)

    def test_radial_aperture(self):
        """Check we cut appropriately with radial aperture"""
        plane_runner = OutputPlaneRunner(100, 0.001, 0.0, "INTERPOLATION", self.verbose)
        plane_runner.output_plane.placement_style = "CENTRE-NORMAL"
        plane_runner.output_plane.output_filename = "" # defaults to "my_plane"
        plane_runner.output_plane.centre = [0.9, 0.1, 0.0]
        plane_runner.output_plane.normal = [1.0, 1.0, 0.0]
        plane_runner.output_plane.radius = 0.1
        width = 0.1
        plane_runner.distribution_str = f"""9
0.0 0.0 0.0 0.0 0.0 0.0
{(width-0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
{(width+0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
{(-width-0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
{(-width+0.001)*0.5**0.5} 0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 {width-0.001} 0.0
0.0 0.0 0.0 0.0 {width+0.001} 0.0
0.0 0.0 0.0 0.0 {-width-0.001} 0.0
0.0 0.0 0.0 0.0 {-width+0.001} 0.0
"""

        plane_runner.execute_fork() # forked process for memory safety
        lines = self.read_lines(plane_runner, "my_plane")
        if self.verbose:
            print("Reference plane test")
            for line in lines:
                print(line[:-1])
        self.assertEqual(len(lines), 5)

    verbose = 0

if __name__ == "__main__":
    unittest.main()