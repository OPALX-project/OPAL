"""
Module to test track run with scaling ffa lattice

This is fast enough to be a unit test; could go in the integration test 
framework but the integration test framework can't do it (yet).
"""

import os
import math
import unittest

import matplotlib
import matplotlib.pyplot

import pyopal.objects.minimal_runner
import pyopal.elements.local_cartesian_offset
import pyopal.elements.scaling_ffa_magnet
import pyopal.elements.asymmetric_enge
import pyopal.objects.field

class ScalingFFARunner(pyopal.objects.minimal_runner.MinimalRunner):
    """Module to test track run (and by extension the entire PyOpal workflow)"""
    def __init__(self):
        """Set up the test"""
        super().__init__()
        self.r0 = 4.0
        self.run_name = "ScalingFFA"
        self.plot_dir = os.getcwd()
        self.verbose = 0

    def build_ffa_magnets(self):
        f_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        lambda_short = 0.07
        lambda_long = 0.14
        spiral_angle = 45
        f_magnet.r0 = self.r0
        f_magnet.b0 = -0.4067
        f_magnet.field_index = 8.0095
        f_magnet.tan_delta = math.tan(math.radians(spiral_angle))
        f_magnet.max_vertical_power = 3
        f_magnet.end_field_model = "f_enge"
        f_magnet.radial_neg_extent = 1.0
        f_magnet.radial_pos_extent = 1.0
        f_magnet.azimuthal_extent = 2.0
        f_magnet.magnet_start = 3*math.pi/160 # note this is a length [m]

        print("ASYMMETRIC ENGE 1")
        f_end = pyopal.elements.asymmetric_enge.AsymmetricEnge()
        print("ASYMMETRIC ENGE 2")
        f_end.set_opal_name("f_enge")
        print("ASYMMETRIC ENGE 2a")
        f_end.lambda_start = 0.14
        f_end.lambda_end = 0.07
        print("ASYMMETRIC ENGE 2b")
        f_end.x0_start = math.pi/20
        f_end.x0_end = f_end.x0_start
        print("ASYMMETRIC ENGE 3")
        f_end.coefficients_start = [0.0, 3.91*cos(math.radians(spiral_angle))]
        f_end.coefficients_end = f_end.coefficients_start
        print("ASYMMETRIC ENGE")
        f_end.update()
        print("ASYMMETRIC ENGE")
        f_magnet.update_end_field()
        print("ASYMMETRIC ENGE")
        """
        d_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        d_magnet.b0 = -0.25592699667
        d_magnet.r0 = f_magnet.r0
        d_magnet.tan_delta = f_magnet.tan_delta
        d_magnet.field_index = f_magnet.field_index
        d_magnet.max_vertical_power = f_magnet.max_vertical_power
        d_magnet.end_field_model = "d_enge"
        d_magnet.radial_neg_extent = f_magnet.radial_neg_extent
        d_magnet.radial_pos_extent = f_magnet.radial_pos_extent

        d_end = pyopal.elements.enge.Enge()
        d_end.set_opal_name("d_enge")
        d_end.x0 = 0.05*d_magnet.r0
        d_end.enge_lambda = 0.15
        d_end.coefficients = f_end.coefficients
        d_end.update()
        d_magnet.update_end_field()
        """
        return [f_magnet] #, d_magnet

    def make_element_iterable(self):
        return self.build_ffa_magnets()

    def field_map_cylindrical(self):
        r_points = [self.r0+i*0.01 for i in range(-100, 100+1, 5)]
        phi_points = [i*0.1 for i in range(-50, 150+1)]
        phi_bins = self.binner(phi_points)
        r_grid = []
        phi_grid = []
        by_grid = []
        for radius in r_points:
            for phi in phi_points:
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
        r_bins = self.binner(r_points)
        phi_bins = self.binner(phi_points)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.hist2d(phi_grid, r_grid, bins=[phi_bins, r_bins], weights=by_grid)
        axes.set_xlabel("r [m]")
        axes.set_ylabel("phi [deg]")
        axes.set_title("by [T]")
        fig_fname = os.path.join(self.plot_dir, "scaling_ffa_map_cyl.png")
        figure.savefig(fig_fname)
        print("Generated cartesian field map in", fig_fname)

    def binner(self, a_list):
        return [a_list[0]+(a_list[1]-a_list[0])*(i-0.5) \
                                                  for i in range(len(a_list)+1)]

    def field_map_cartesian(self):
        x_points = [i*0.03 for i in range(-0, 200+1)]
        y_points = [i*0.03 for i in range(-0, 200+1)]
        x_grid = []
        y_grid = []
        by_grid = []
        for x in x_points:
            for y in y_points:
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
        x_bins = self.binner(x_points)
        y_bins = self.binner(y_points)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        hist = axes.hist2d(x_grid, y_grid, bins=[x_bins, y_bins], weights=by_grid)
        axes.set_xlabel("x [m]")
        axes.set_ylabel("y [m]")
        axes.set_title("by [T]")
        figure.colorbar(hist[3])
        fig_fname = os.path.join(self.plot_dir, "scaling_ffa_map_cart.png")
        figure.savefig(fig_fname)
        print("Generated cartesian field map in", fig_fname)

    def field_maps(self):
        self.field_map_cartesian()
        self.field_map_cylindrical()


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
        if self.do_verbose:
            runner.verbose = 3
            runner.postprocess = runner.field_maps
        runner.run_one_fork()

    do_verbose = False

def main():
    runner = ScalingFFARunner()
    runner.verbose = 3
    runner.postprocess = runner.field_maps
    runner.run_one_fork()


if __name__ == "__main__":
    main()
#    TestTrackRun.do_verbose = True # verbose if called from command line
#    unittest.main()
