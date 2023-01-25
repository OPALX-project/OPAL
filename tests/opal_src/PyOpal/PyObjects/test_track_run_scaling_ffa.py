"""
Module to test track run with scaling ffa lattice

This is fast enough to be a unit test; could go in the integration test 
framework but the integration test framework can't do it (yet).
"""

import os
import math
import unittest

import pyopal.objects.minimal_runner
import pyopal.elements.local_cartesian_offset
import pyopal.elements.scaling_ffa_magnet
import pyopal.elements.asymmetric_enge
import pyopal.elements.probe
import pyopal.objects.field

C_LIGHT = 299792458 # m/s

class ScalingFFARunner(pyopal.objects.minimal_runner.MinimalRunner):
    """Module to test simulation of a scaling ffa (and by extension the entire PyOpal workflow)"""
    #UPDATE OFFSETS SO THAT THEY WORK!

    def __init__(self):
        """Set up the test"""
        super().__init__()
        self.r0 = 4.0
        self.run_name = "ScalingFFA"
        self.plot_dir = os.getcwd()
        self.verbose = 0
        self.probe_id = 1
        self.n_cells = 16
        self.bend_direction = -1 # set to -1 for ring centred on x,y = 8,0 m
        ke = 3e-3 # 3 MeV
        self.momentum = ((ke+self.mass)**2-self.mass**2)**0.5
        self.cell_length = 2*math.pi/self.n_cells*self.r0 # [m]
        self.max_steps = 1000
        self.set_number_of_steps(1e-2)

    def set_number_of_steps(self, step_size):
        global C_LIGHT
        beta_relativistic = self.momentum/(self.momentum**2+self.mass**2)**0.5
        time_step = step_size/(beta_relativistic*C_LIGHT) # seconds
        self.steps_per_turn = 1000
        self.time_per_turn = time_step*1000 # [seconds], in fact this is the time to do 1000 steps
        print("beta_relativistic", beta_relativistic, self.time_per_turn, 1/self.time_per_turn)

    def build_ffa_magnets(self):
        f_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        lambda_short = 0.07
        lambda_long = 0.14
        spiral_angle = 45
        f_magnet.r0 = self.bend_direction*self.r0
        f_magnet.b0 = -0.4067
        f_magnet.field_index = 8.0095
        f_magnet.tan_delta = math.tan(math.radians(spiral_angle))
        f_magnet.max_vertical_power = 3
        f_magnet.end_field_model = "f_enge"
        f_magnet.radial_neg_extent = 1.0
        f_magnet.radial_pos_extent = 1.0
        # azimuthal extent [m] defines the extent of the bounding box
        f_magnet.azimuthal_extent = self.cell_length
        # magnet start [m] defines where the fringe field starts rising
        # relative to the element start
        f_magnet.magnet_start = self.cell_length/4
        # magnet end [m] defines where this element ends i.e. placement of the
        # next element
        f_magnet.magnet_end = self.cell_length/2

        # f_end defines the actual value of the field as a function of phi
        f_end = pyopal.elements.asymmetric_enge.AsymmetricEnge()
        f_end.set_opal_name("f_enge")
        f_end.lambda_start = lambda_long
        f_end.lambda_end = lambda_short
        f_end.x0_start = math.pi/20
        f_end.x0_end = f_end.x0_start
        f_end.coefficients_start = [0.0, 3.91*math.cos(math.radians(spiral_angle))]
        f_end.coefficients_end = f_end.coefficients_start
        f_end.update()
        f_magnet.update_end_field()

        d_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        d_magnet.b0 = 0.25592699667
        d_magnet.r0 = f_magnet.r0
        d_magnet.tan_delta = f_magnet.tan_delta
        d_magnet.field_index = f_magnet.field_index
        d_magnet.max_vertical_power = f_magnet.max_vertical_power
        d_magnet.end_field_model = "d_enge"
        d_magnet.radial_neg_extent = f_magnet.radial_neg_extent
        d_magnet.radial_pos_extent = f_magnet.radial_pos_extent
        d_magnet.azimuthal_extent = self.cell_length
        d_magnet.magnet_start = self.cell_length/8 # note this is a length [m]
        d_magnet.magnet_end = self.cell_length/4 # note this is a length [m]

        d_end = pyopal.elements.asymmetric_enge.AsymmetricEnge()
        d_end.set_opal_name("d_enge")
        d_end.x0_start = math.pi/40
        d_end.x0_end = math.pi/40
        d_end.lambda_start = lambda_short
        d_end.lambda_end = lambda_short
        d_end.coefficients_start = f_end.coefficients_start
        d_end.coefficients_end = f_end.coefficients_end
        d_end.update()
        d_magnet.update_end_field()

        return [f_magnet, d_magnet]

    def build_probe(self, phi_degrees):
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
        # delta_angle is the total length taken up by the f_magnet and the 
        # d_magnet together, along a line having r = self.r0
        delta_length = (f_magnet.magnet_end+d_magnet.magnet_end) # [m]
        # required angle is the angle that the drift is required to fill
        required_angle = (self.cell_length-delta_length)/self.r0
        drift = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        drift.end_position_x=self.r0*math.sin(required_angle)
        drift.end_position_y=-self.bend_direction*self.r0*(1-math.cos(required_angle))
        drift.end_normal_x=math.cos(required_angle)
        drift.end_normal_y=-self.bend_direction*math.sin(required_angle)
        return drift


    def make_element_iterable(self):
        probes = [self.build_probe(360.0/self.n_cells*i) for i in range(self.n_cells)]
        f_magnet, d_magnet = self.build_ffa_magnets()
        drift = self.build_drift_space(f_magnet, d_magnet)
        cell = [f_magnet, d_magnet, drift]*16
        return cell+probes

    def plots(self, plot_dir = None):
        mapper = FieldMapper()
        mapper.plot_dir = self.plot_dir
        mapper.load_tracks(os.path.join(self.tmp_dir, self.run_name+"-trackOrbit.dat"))

        mapper.r_points = [self.r0+i*0.01 for i in range(-120, 120+1, 5)]
        mapper.phi_points = [i*0.1 for i in range(0, 1800+1)]
        mapper.field_map_cylindrical()

        mapper.x_points = [i*0.03 for i in range(-200, 800+1)]
        mapper.y_points = [i*0.03 for i in range(-800, 800+1)]
        mapper.field_map_cartesian()


class FieldMapper(object):
    def __init__(self):
        self.r_points = []
        self.phi_points = []

        self.x_points = []
        self.y_points = []

        self.verbose = 0
        self.cmap = "PiYG"
        self.plot_dir = os.getcwd()

    @classmethod
    def binner(self, a_list):
        return [a_list[0]+(a_list[1]-a_list[0])*(i-0.5) \
                                                  for i in range(len(a_list)+1)]

    def load_tracks(self, track_orbit):
        with open(track_orbit) as fin:
            for line in range(2):
                print(fin.readline())

    def gen_cmap(self, by_grid):
        min_by = min(by_grid)
        max_by = max(by_grid)
        cmax = max(abs(min_by), abs(max_by))
        return min_by, max_by, cmax

    def field_map_cylindrical(self):
        r_grid = []
        phi_grid = []
        by_grid = []
        for radius in self.r_points:
            for phi in self.phi_points:
                r_grid.append(radius)
                phi_grid.append(phi)
                point = (radius*math.cos(math.radians(phi)),
                         radius*math.sin(math.radians(phi)),
                         0,
                         0)
                value = pyopal.objects.field.get_field_value(*point)
                by_grid.append(value[3])
                if self.verbose > 0:
                    print("Field value at r, phi", radius, round(phi, 2),
                          "point", point,
                          "is B:", value[1:4],
                          "E:", value[4:])
        try:
            import matplotlib
            import matplotlib.pyplot
        except ImportError:
            print("Matplotlib not imported - not making plots")
            return
        r_bins = self.binner(self.r_points)
        phi_bins = self.binner(self.phi_points)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        min_by, max_by, cmax = self.gen_cmap(by_grid)
        axes.hist2d(phi_grid, r_grid, bins=[phi_bins, r_bins], weights=by_grid,
                    cmin=min_by, cmax=max_by, cmap=self.cmap, vmin=-cmax, vmax=cmax)
        axes.set_xlabel("r [m]")
        axes.set_ylabel("phi [deg]")
        axes.set_title("by [T]")
        fig_fname = os.path.join(self.plot_dir, "scaling_ffa_map_cyl.png")
        figure.savefig(fig_fname)
        print("Generated cylindrical field map in", fig_fname)
        return figure

    def field_map_cartesian(self):
        x_grid = []
        y_grid = []
        by_grid = []
        for x in self.x_points:
            for y in self.y_points:
                x_grid.append(x)
                y_grid.append(y)
                point = (x, y, 0, 0)
                value = pyopal.objects.field.get_field_value(*point)
                by_grid.append(value[3])
                if self.verbose > 0:
                    print("Field value at point", point,
                          "is B:", value[1:4], "E:", value[4:])
        try:
            import matplotlib
            import matplotlib.pyplot
        except ImportError:
            print("Matplotlib not imported - not making plots")
            return
        x_bins = self.binner(self.x_points)
        y_bins = self.binner(self.y_points)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        min_by, max_by, cmax = self.gen_cmap(by_grid)
        hist = axes.hist2d(x_grid, y_grid, bins=[x_bins, y_bins], weights=by_grid,
                    cmin=min_by, cmax=max_by, cmap=self.cmap, vmin=-cmax, vmax=cmax)
        axes.set_xlabel("x [m]")
        axes.set_ylabel("y [m]")
        axes.set_title("by [T]")
        figure.colorbar(hist[3])
        fig_fname = os.path.join(self.plot_dir, "scaling_ffa_map_cart.png")
        figure.savefig(fig_fname)
        print("Generated cartesian field map in", fig_fname)
        return figure


    distribution_str = """1
3.944586177309523 -0.02776333011661966 0.0 -0.0049890385556281445 0.1584654928597547 -0.0016918209895814252
"""

class TestTrackRun(unittest.TestCase):
    """Test class for track_run"""
    def test_run_one(self):
        """
        Test that we can run okay without an exception.

        If running from the command line, will spit out the OPAL log to screen
        """
        runner = ScalingFFARunner()
        runner.plot_dir = os.getcwd()
        if self.do_verbose:
            runner.verbose = 0
            runner.postprocess = runner.plots
        runner.execute_fork()
        print("Finished\n\n")

    do_verbose = False

if __name__ == "__main__":
    TestTrackRun.do_verbose = True # verbose if called from command line
    unittest.main()
