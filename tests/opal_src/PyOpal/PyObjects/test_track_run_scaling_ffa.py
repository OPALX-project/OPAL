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
import pyopal.elements.local_cartesian_offset
import pyopal.elements.scaling_ffa_magnet
import pyopal.elements.asymmetric_enge
import pyopal.elements.probe
import pyopal.objects.field
import pyopal.objects.minimal_runner
try:
    # for plotting, we require matplotlib. If it cannot be imported, plotting
    # will not be performed.
    import matplotlib
    matplotlib.use('agg') # disable matplot gui during this test
    import pyopal.objects.ffa_field_mapper
    HAS_FFA_MAPPER = True
except ImportError:
    HAS_FFA_MAPPER = False

C_LIGHT = 299792458 # m/s

class ScalingFFARunner(pyopal.objects.minimal_runner.MinimalRunner):
    """
    Module to test simulation of a scaling ffa (and by extension the entire 
    PyOpal workflow)
    """
    def __init__(self):
        """Here we set some basic lattice parameters"""
        super().__init__()
        self.r0 = 4.0
        self.run_name = "ScalingFFA"
        self.plot_dir = os.getcwd()
        self.verbose = 0
        self.probe_id = 1
        self.n_cells = 16
        self.bend_direction = 1 # set to -1 for ring centred on x,y = 8,0 m
        self.dr = 0.5
        self.spiral_angle = math.pi/6
        ke = 3e-3 # 3 MeV
        self.momentum = ((ke+self.mass)**2-self.mass**2)**0.5
        self.cell_length = 2*math.pi/self.n_cells*self.r0 # [m]
        self.f_start = self.cell_length/8.0
        self.f_end_length = self.cell_length/32.0
        self.f_centre_length = self.cell_length/4.0
        self.f_end = (self.f_start+self.f_centre_length+self.f_end_length*4)
        self.max_steps = 10
        self.set_number_of_steps(1e-2)

    def set_number_of_steps(self, step_size):
        """
        Do some maths to set the number of steps using the step size as input
        - step_size: step size in metres
        """
        global C_LIGHT
        beta_relativistic = self.momentum/(self.momentum**2+self.mass**2)**0.5
        time_step = step_size/(beta_relativistic*C_LIGHT) # seconds
        self.steps_per_turn = 1000
        self.time_per_turn = time_step*1000 # [seconds], in fact this is the time to do 1000 steps
        print("beta_relativistic", beta_relativistic, self.time_per_turn, 1/self.time_per_turn)

    def build_ffa_magnets(self):
        """
        The lattice has two FFA magnets, a (F)ocusing magnet and a (D)efocusing
        magnet. In this example, the f magnet uses the default end field model,
        which is a tanh model. The d magnet uses a more general enge field
        model.
        """
        # f magnet is an example of a magnet using the default tanh end field
        f_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        f_magnet.set_opal_name("f_magnet")
        f_magnet.r0 = self.bend_direction*self.r0
        f_magnet.b0 = -0.5
        f_magnet.field_index = 8.0095
        f_magnet.tan_delta = math.tan(self.spiral_angle)
        f_magnet.max_vertical_power = 3
        f_magnet.radial_neg_extent = self.dr/2
        f_magnet.radial_pos_extent = self.dr/2
        # azimuthal extent [m] defines the extent of the bounding box
        f_magnet.azimuthal_extent = self.cell_length
        # magnet start [m] defines where the fringe field starts rising
        # relative to the element start
        f_magnet.magnet_start = self.f_start
        # magnet end length [m] defines how quickly the fringe field drops
        f_magnet.end_length = self.f_end_length
        # magnet centre length [m] defines the length of the flat top
        f_magnet.centre_length = self.f_centre_length
        # magnet end [m] defines where the next element will be placed
        f_magnet.magnet_end = self.f_end

        # d magnet is an example of a magnet using an end field "plug in" module
        d_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        d_magnet.set_opal_name("d_magnet")
        d_magnet.b0 = 0.25
        d_magnet.r0 = f_magnet.r0
        d_magnet.tan_delta = f_magnet.tan_delta
        d_magnet.field_index = f_magnet.field_index
        d_magnet.max_vertical_power = f_magnet.max_vertical_power
        d_magnet.end_field_model = "d_enge"
        d_magnet.radial_neg_extent = f_magnet.radial_neg_extent
        d_magnet.radial_pos_extent = f_magnet.radial_pos_extent
        d_magnet.azimuthal_extent = self.cell_length
        d_magnet.magnet_start = f_magnet.magnet_start
        d_magnet.magnet_end = f_magnet.magnet_end/2

        d_end = pyopal.elements.asymmetric_enge.AsymmetricEnge()
        d_end.set_opal_name("d_enge")
        d_end.x0_start = self.f_centre_length/4
        d_end.x0_end = self.f_centre_length/4
        d_end.lambda_start = self.f_end_length/2
        d_end.lambda_end = self.f_end_length/2
        d_end.coefficients_start = [0.0, 1.0]
        d_end.coefficients_end = [0.0, 1.0]
        d_end.update()
        d_magnet.update_end_field()
        print(f"Built magnets with r0 {f_magnet.r0}")

        return [f_magnet, d_magnet]

    def build_probe(self, phi_degrees):
        """
        Build a single probe, positioned along a radius at a given azimuthal 
        angle
        - phi_degrees: the azimuthal angle at which the probe is placed.
        """
        name = "ring_probe_"+str(self.probe_id).rjust(3, "0")
        self.probe_id += 1
        probe = pyopal.elements.probe.Probe()
        probe.x_start = 0.0
        probe.y_start = 0.0
        probe.x_end = self.r0*2*1000*math.cos(math.radians(phi_degrees))
        probe.y_end = self.r0*2*1000*math.sin(math.radians(phi_degrees))
        probe.output_filename = name
        probe.set_opal_name(name)
        return probe

    def build_drift_space(self, f_magnet, d_magnet):
        """
        Build a drift space at the end of the FD magnet pair to make the cell up
        to self.cell_length, in the local coordinate system
        """
        # delta_angle is the total length taken up by the f_magnet and the 
        # d_magnet together, along a line having r = self.r0
        delta_length = (f_magnet.magnet_end+d_magnet.magnet_end) # [m]
        # required angle is the angle that the drift is required to fill
        required_angle = (self.cell_length-delta_length)/self.r0
        drift = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        drift.end_position_x=self.bend_direction*self.r0*(math.cos(required_angle)-1)
        drift.end_position_y=self.r0*math.sin(required_angle)
        drift.end_normal_x=-self.bend_direction*math.sin(required_angle)
        drift.end_normal_y=math.cos(required_angle)
        return drift


    def make_element_iterable(self):
        """
        Overload the method in MinimalRunner to place the field elements in the
        lattice.
        """
        probes = [self.build_probe(360.0/self.n_cells*i) for i in range(self.n_cells)]
        f_magnet, d_magnet = self.build_ffa_magnets()
        drift = self.build_drift_space(f_magnet, d_magnet)
        #cell = [f_magnet]*3 #, drift]*3
        cell = [f_magnet, d_magnet, drift]*16
        return cell+probes

    def plots(self, plot_dir = None):
        """
        Make plots showing rectangular and cylindrical field maps. The
        granularity and region mapped can be adjusted by changing r_points,
        phi_points (cylindrical) and x_point, y_points (cartesian)
        """
        global HAS_FFA_MAPPER
        if not HAS_FFA_MAPPER: # matplotlib dependency
            return
        mapper = pyopal.objects.ffa_field_mapper.FFAFieldMapper()
        mapper.plot_dir = self.plot_dir
        mapper.load_tracks(os.path.join(self.tmp_dir, self.run_name+"-trackOrbit.dat"))

        mapper.radial_contours = [
            {"radius":self.r0, "linestyle":"dashed", "colour":"grey", "label":"r0"},
            {"radius":self.r0-self.dr/2, "linestyle":"-", "colour":"grey", "label":"r0-\nradial_neg_extent"},
            {"radius":self.r0+self.dr/2, "linestyle":"-", "colour":"grey", "label":"r0+\nradial_pos_extent"},
        ]


        phi_start = math.degrees(self.f_start/self.r0)
        phi_middle = math.degrees((self.f_start+self.f_centre_length/2)/self.r0)
        phi_fringe_end = math.degrees((self.f_start+self.f_centre_length)/self.r0)
        phi_end = math.degrees(self.f_end/self.r0)
        spiral_deg = math.degrees(self.spiral_angle)
        mapper.spiral_contours = [
            {"phi0":phi_start, "r0":self.r0, "spiral_angle":spiral_deg, "linestyle":"dashed", "colour":"blue", "label":"magnet_start/r0"},
            {"phi0":phi_middle, "r0":self.r0, "spiral_angle":spiral_deg, "linestyle":"dashed", "colour":"grey", "label":"(magnet_start+\ncentre_length/2)/r0"},
            {"phi0":phi_fringe_end, "r0":self.r0, "spiral_angle":spiral_deg, "linestyle":"dashed", "colour":"grey", "label":"(magnet_start+\ncentre_length)/r0"},
            {"phi0":phi_end, "r0":self.r0, "spiral_angle":spiral_deg, "linestyle":"dashed", "colour":"blue", "label":"magnet_end/r0"},
        ]

        mapper.r_points = [self.r0+i*0.001 for i in range(-1200, 1200+1, 10)]
        mapper.phi_points = [i/8.0/10 for i in range(0, 1800+1, 5)]
        mapper.field_map_cylindrical()

        mapper.x_points = [i*0.05 for i in range(-100, 100+1, 5)]
        mapper.y_points = [i*0.05 for i in range(-100, 100+1, 5)]

        mapper.field_map_cartesian()
        mapper.oned_field_map(self.r0)


        n_elements = pyopal.objects.field.get_number_of_elements()
        if not self.verbose:
            return
        print("PyOpal built", n_elements, "elements, not including the RingDefinition:")
        for i in range(pyopal.objects.field.get_number_of_elements()):
            start_pos = pyopal.objects.field.get_element_start_position(i)
            start_rad = (start_pos[0]**2+start_pos[1]**2)**0.5
            start_phi = math.atan2(start_pos[1], start_pos[0])
            print("   ", pyopal.objects.field.get_element_name(i),
                  "starting at radius", start_rad, "m and angle",
                  start_phi, "rad")

    distribution_str = """1
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
"""

class TestTrackRun(unittest.TestCase):
    """Test class for track_run"""
    def test_run_one(self):
        """
        Test that we can run okay without an exception.

        If running from the command line, will spit out the OPAL log to screen.
        """
        runner = ScalingFFARunner()
        runner.bend_direction = 1
        runner.plot_dir = os.getcwd()
        runner.postprocess = runner.plots
        runner.execute_fork()
        print("Finished\n\n")

    def test_end_length(self):
        """
        Test that we can run okay without an exception.

        If running from the command line, will spit out the OPAL log to screen.
        """
        runner = ScalingFFARunner()
        runner.bend_direction = 1
        runner.plot_dir = os.getcwd()
        runner.postprocess = runner.plots
        runner.execute_fork()
        print("Finished\n\n")



if __name__ == "__main__":
    unittest.main()
