"""
Minimal set of methods to build an OPAL simulation and run it.
This takes about 0.1 s to run on my average desktop PC.
"""
import os
import tempfile

import pyopal.objects.track_run
import pyopal.objects.beam
import pyopal.objects.distribution
import pyopal.objects.line
import pyopal.elements.ring_definition
import pyopal.elements.local_cartesian_offset
import pyopal.objects.field_solver
import pyopal.objects.track
import pyopal.objects.option

class MinimalRunner(object):
    """Class to run a minimal OPAL setup.

    Individual methods make_<BLAH> handle set up of each OPAL "global" object:
    - self.field_solver
    - self.distribution
    - self.beam
    - self.track
    - self.track_run
    - self.line
    - self.ring
    """
    def __init__(self):
        """Initialise to empty data"""
        self.field_solver = None
        self.distribution = None
        self.line = None
        self.ring = None
        self.beam = None
        self.track = None
        self.run = None
        self.track_run = None
        self.option = None

        self.tmp_dir = tempfile.mkdtemp()
        self.distribution_filename = os.path.join(self.tmp_dir,
                                                  "distribution.dat")
        self.r0 = 1.0
        self.momentum = 0.1
        self.mass = 0.93827208816
        self.run_name = None

    def make_field_solver(self):
        """Make an empty fieldsolver

        The fieldsolver has the job of solving the space charge problem on
        successive time steps. In this example, the FieldSolver is switched off
        (i.e. set to type = "NONE").
        """
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

    def make_distribution(self):
        """Make a distribution

        The distribution is the initial beam distribution of the PyOPAL
        simulation. In this example, a distribution loaded from a tempfile is
        called, where the tempfile is written to disk dynamically at runtime.
        """
        dist_file = open(self.distribution_filename, "w+")
        dist_file.write(self.distribution_str)
        dist_file.flush()
        dist_file.close()
        self.distribution = pyopal.objects.distribution.Distribution()
        self.distribution.set_opal_name("SuperDist")
        self.distribution.type = "FROMFILE"
        self.distribution.filename = self.distribution_filename
        self.distribution.register()
        return self.distribution

    def make_beam(self):
        """Make a beam

        The beam holds the global/default beam distribution information, such as
        the mass of particles in the beam, number of particles in the beam and
        so on.
        """
        beam = pyopal.objects.beam.Beam()
        beam.set_opal_name("SuperBeam")
        beam.mass = self.mass
        beam.momentum = self.momentum
        beam.charge = 1.0
        beam.beam_frequency = 1.0
        beam.number_of_slices = 10
        beam.number_of_particles = 1
        beam.register()
        self.beam = beam

    @classmethod
    def null_drift(cls):
        """Returns a drift of length 0

        OPAL requires at least one element in each beam line. For this simplest
        example a drift of length 0 is used.
        """
        drift = pyopal.elements.local_cartesian_offset.LocalCartesianOffset()
        drift.end_position_x=0.0
        drift.end_position_y=0.0
        drift.end_normal_x=1.0
        drift.end_normal_y=0.0
        return drift

    def make_ring(self):
        """Make a RingDefinition object.

        The RingDefinition holds default parameters for a ring initial 
        conditions, in particular the initial cylindrical coordinates for the 
        first element placement and beam, and the minimum and maximum radius
        allowed before particles are considered lost. The ring can be appended
        to self.line and used with OPAL cyclotron mode.
        """
        self.ring = pyopal.elements.ring_definition.RingDefinition()
        self.ring.set_opal_name("a_ring")
        self.ring.lattice_initial_r = self.r0
        self.ring.beam_initial_r = self.r0
        self.ring.minimum_r = self.r0/2
        self.ring.maximum_r = self.r0*2
        self.ring.is_closed = False

    def make_option(self):
        """Options enable setting of global control variables.

        For a full list of variables see the Option docs.
        """
        self.option = pyopal.objects.option.Option()
        self.option.echo = False
        self.option.info = False
        self.option.spt_dump_frequency = 2.0
        self.option.execute()

    def make_element_iterable(self):
        """
        Return an iterable (e.g. list) of elements to append to the line

        By default, returns an empty list. User can overload this method.
        """
        return []

    def make_line(self):
        """Make a Line object.

        The Line holds a sequence of beam elements.
        """
        self.line = pyopal.objects.line.Line()
        self.line.set_opal_name("test_line")
        try:
            self.line.append(self.ring)
        except Exception:
            print(self.ring_error)
            raise
        self.line.append(self.null_drift())
        an_element_iter = self.make_element_iterable()
        for element in an_element_iter:
            self.line.append(element)
        self.line.register()   

    def make_track(self):
        """Make a track object.

        The track object handles the interface between tracking and the beam
        elements.
        """
        track = pyopal.objects.track.Track()
        track.line = "test_line"
        track.beam = "SuperBeam"
        self.track = track
        self.track.execute()

    def make_track_run(self):
        """Make a TrackRun

        The TrackRun handles the interface between the Track, distribution, 
        field solver and calls the actual tracking routines.
        """
        run = pyopal.objects.track_run.TrackRun()
        run.method = "CYCLOTRON-T"
        run.keep_alive = True
        run.beam_name = "SuperBeam"
        run.distribution = ["SuperDist"]
        run.field_solver = "FIELDSOLVER"
        run.steps_per_turn = 100
        self.track_run = run

    def run_one(self):
        """Set up and run a simulation"""
        here = os.getcwd()
        try:
            os.chdir(self.tmp_dir)
            self.make_option()
            self.make_distribution()
            self.make_field_solver()
            self.make_beam()
            self.make_ring()
            self.make_line()
            self.make_track()
            self.make_track_run()
            if self.run_name:
                self.track_run.set_run_name(self.run_name)
            self.track_run.execute()
        except:
            raise
        finally:
            print("Finished running in directory", os.getcwd())
            os.chdir(here)

    distribution_str = """1
0.0 0.0 0.0 0.0 0.0 0.0
"""

    ring_error = "Failed to append ring to the line. Try running make_ring()"+\
    " before calling make_line()."

def main():
    opal = MinimalRunner()
    opal.run_one()

if __name__ == "__main__":
    main()
