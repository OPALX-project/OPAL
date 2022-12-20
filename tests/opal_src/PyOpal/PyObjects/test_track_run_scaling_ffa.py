"""
Module to test track run with scaling ffa lattice

This is fast enough to be a unit test; could go in the integration test 
framework but the integration test framework can't do it yet.
"""

import math
import os
import sys
import tempfile
import subprocess
import unittest

import pyopal.elements.scaling_ffa_magnet
import pyopal.objects.track_run
import pyopal.objects.beam
import pyopal.objects.distribution
import pyopal.objects.line
import pyopal.elements.ring_definition
import pyopal.elements.local_cartesian_offset
import pyopal.objects.field_solver
import pyopal.objects.track
import pyopal.objects.parser
import pyopal.elements.enge
import pyopal.objects.field


class TrackRunExecute():
    """Module to test track run (and by extension the entire PyOpal workflow)"""
    def __init__(self):
        """Set up the test"""
        self.here = os.getcwd()
        self.tmp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.tmp_dir.name)
        self.field_solver = None
        self.line = None
        self.ring = None
        self.offset = None
        self.distribution = None
        self.distribution_file = tempfile.NamedTemporaryFile("w+")
        self.track_run = pyopal.objects.track_run.TrackRun()
        self.r0 = 4.0


    def make_field_solver(self):
        """Make an empty fieldsolver"""
        self.field_solver = pyopal.objects.field_solver.FieldSolver()
        self.field_solver.type = "NONE"
        self.field_solver.mesh_size_x = 5
        self.field_solver.mesh_size_y = 5
        self.field_solver.mesh_size_t = 5
        self.field_solver.parallelize_x = False
        self.field_solver.parallelize_y = False
        self.field_solver.parallelize_t = False
        self.field_solver.boundary_x = "open"
        self.field_solver.boundary_y = "open"
        self.field_solver.boundary_t = "open"
        self.field_solver.bounding_box_increase = 2
        self.field_solver.register()

    @classmethod
    def make_drift(cls):
        """Returns a drift of length 0"""
        drift = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        drift.end_position_x=0.0
        drift.end_position_y=0.0
        drift.end_normal_x=1.0
        drift.end_normal_y=0.0
        return drift

    def make_line(self):
        """Make a beamline, just consisting of a drift section"""
        drift = self.make_drift()
        self.line = pyopal.objects.line.Line()
        self.ring = pyopal.elements.ring_definition.RingDefinition()
        self.ring.set_opal_name("a_ring")
        self.line.set_opal_name("ffa_line")
        self.ring.lattice_initial_r = self.r0
        self.ring.beam_initial_r = 0.0
        self.ring.minimum_r = 0.5
        self.ring.maximum_r = 10.0
        self.ring.is_closed = False
        self.offset = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        self.offset.end_position_x = 0.0
        self.offset.end_position_y = 1.0
        self.offset.normal_x = 1.0

        self.line.append(self.ring)
        f_magnet, d_magnet = self.make_scaling_ffa_magnets()
        self.line.append(drift)
        self.line.append(f_magnet)
        #self.line.append(d_magnet)
        self.line.register()

    def make_distribution(self):
        """Make a distribution, from tempfile data"""
        self.distribution_file.write(self.distribution_str)
        self.distribution_file.flush()
        self.distribution = pyopal.objects.distribution.Distribution()
        self.distribution.set_opal_name("SuperDist")
        self.distribution.type = "FROMFILE"
        self.distribution.fname = self.distribution_file.name
        self.distribution.register()

        return self.distribution

    def run_one(self):
        """Set up and run a simulation"""
        self.make_line()
        self.make_distribution()
        self.make_field_solver()

        beam = pyopal.objects.beam.Beam()
        beam.set_opal_name("SuperBeam")
        beam.mass = 0.93827208816
        beam.charge = 1.0
        beam.momentum = 0.1
        beam.beam_frequency = 1.0
        beam.number_of_slices = 10
        beam.number_of_particles = 1
        beam.register()

        track = pyopal.objects.track.Track()
        track.line = "ffa_line"
        track.beam = "SuperBeam"
        run = pyopal.objects.track_run.TrackRun()
        run.method = "CYCLOTRON-T"
        run.keep_alive = True
        run.beam_name = "SuperBeam"
        run.distribution = ["SuperDist"]
        run.field_solver = "FIELDSOLVER"
        track.execute()
        run.execute()
        track.execute()
        run.execute()
        self.field_map_cartesian()

    def field_map_cylindrical(self):
        r_points = [self.r0+i*0.01 for i in range(1, 1+1)]
        phi_points = [i*1 for i in range(-5, 15+1)]
        phi_bins = [phi_points[0]+(phi_points[1]-phi_points[0])*(i-0.5) for i in range(len(phi_points)+1)]
        r_grid = []
        phi_grid = []
        by_grid = []
        for radius in r_points:
            for phi in phi_points:
                r_grid.append(radius)
                phi_grid.append(phi)
                point = (radius*math.cos(math.radians(phi)), radius*math.sin(math.radians(phi)), 0, 0)
                value = pyopal.objects.field.get_field_value(*point)
                by_grid.append(value[3])
                print("Field value at r, phi", radius, round(phi, 2), "point", point, "is B:", value[1:4], "E:", value[4:])
        if TestTrackRun.verbose < 3:
            return
        try:
            import matplotlib
            import matplotlib.pyplot
        except ImportError:
            print("Matplotlib not imported - not making plots")
            return
        binner = lambda a_list: [a_list[0]+(a_list[1]-a_list[0])*(i-0.5) for i in range(len(a_list)+1)]
        r_bins = binner(r_points)
        phi_bins = binner(phi_points)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.hist2d(r_grid, phi_grid, bins=[r_bins, phi_bins], weights=by_grid)
        axes.set_xlabel("r [m]")
        axes.set_ylabel("phi [deg]")
        axes.set_title("by [T]")
        figure.savefig(os.path.expandvars("${HOME}/fig.png"))

    def make_scaling_ffa_magnets(self):
        f_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        f_magnet.r0 = self.r0
        f_magnet.b0 = 0.5
        f_magnet.field_index = 1
        f_magnet.tan_delta = math.tan(math.radians(10.0))
        f_magnet.max_vertical_power = 3
        f_magnet.end_field_model = "f_enge"
        f_magnet.radial_neg_extent = 1.0
        f_magnet.radial_pos_extent = 1.0
        f_magnet.azimuthal_extent = 2.0
        f_magnet.magnet_start = 1.0

        f_end = pyopal.elements.enge.Enge()
        f_end.set_opal_name("f_enge")
        f_end.x0 = 0.1
        f_end.enge_lambda = 0.05
        f_end.coefficients = [-1.3, 3.42, -1.2, 0.3]
        f_end.update()
        f_magnet.update_end_field()

        d_magnet = pyopal.elements.scaling_ffa_magnet.ScalingFFAMagnet()
        d_magnet.b0 = -0.25592699667
        d_magnet.r0 = f_magnet.r0
        d_magnet.tan_delta = f_magnet.tan_delta
        d_magnet.field_index = f_magnet.field_index
        d_magnet.max_vertical_power = f_magnet.max_vertical_power
        d_magnet.end_field_model = "d_enge"

        d_end = pyopal.elements.enge.Enge()
        d_end.set_opal_name("d_enge")
        d_end.x0 = 0.05*d_magnet.r0
        d_end.enge_lambda = 0.15
        d_end.coefficients = f_end.coefficients
        d_end.update()
        d_magnet.update_end_field()

        return f_magnet, d_magnet

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
                print("Field value at point", point, "is B:", value[1:4], "E:", value[4:])
        if TestTrackRun.verbose < 3:
            return
        try:
            import matplotlib
            import matplotlib.pyplot
        except ImportError:
            print("Matplotlib not imported - not making plots")
            return
        binner = lambda a_list: [a_list[0]+(a_list[1]-a_list[0])*(i-0.5) for i in range(len(a_list)+1)]
        x_bins = binner(x_points)
        y_bins = binner(y_points)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        hist = axes.hist2d(x_grid, y_grid, bins=[x_bins, y_bins], weights=by_grid)
        axes.set_xlabel("x [m]")
        axes.set_ylabel("y [m]")
        axes.set_title("by [T]")
        figure.colorbar(hist[3])
        figure.savefig(os.path.expandvars("${HOME}/fig.png"))

    def __del__(self):
        """move back to the old cwd"""
        os.chdir(self.here)
        self.tmp_dir.cleanup()
        self.distribution_file.close()


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
        log = tempfile.TemporaryFile("w+")
        proc = subprocess.run(["python", __file__, "run_test_track_run"],
            stdout=log, stderr=subprocess.STDOUT, check=False)
        log.seek(0)
        error_message=""
        if self.verbose > 0:
            error_message = log.read()
            if self.verbose > 1:
                print(error_message)
        log.close()
        self.assertEqual(proc.returncode, 0, msg=error_message)

    verbose = 3

if __name__ == "__main__":
    if "run_test_track_run" in sys.argv:
        TrackRunExecute().run_one()
    else:
        # verbose = 0 -> silent
        # verbose = 1 -> print on error
        # verbose = 2 -> always print
        # verbose = 3 -> make a plot in ${HOME}
        unittest.main()
