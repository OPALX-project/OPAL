# Basic lattice set up defaults - either use as an example or inherit and 
# overload required methods
#
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
#

"""
Minimal set of methods to build an OPAL simulation and run it.
This takes about 0.1 s to run on my average desktop PC.
"""
import os
import sys
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

    Individual methods make_foo handle set up of each OPAL "global" object:
    - self.field_solver
    - self.distribution
    - self.beam
    - self.track
    - self.track_run
    - self.line
    - self.ring
    Explicitly, the OPAL routines are only called in make_foo and run_one so
    that there is no OPAL things initialised and run_one_fork can be safely
    called.

    There are three hooks for user to overload and do stuff:
    - make_element_iterable: add extra elements to the line
    - preprocess: add some stuff to do just before tracking starts
    - postprocess: add some stuff to do just after tracking ends
    user can do postprocessing after calling run_one, but run_one_fork isolates
    memory allocation into a forked process so parent process does not have
    access to memory for e.g. checking details of the field maps, etc.
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
        self.max_steps = 100
        self.r0 = 1.0 # [m]
        self.momentum = 0.1 # [GeV/c]
        self.mass = 0.93827208816 # [GeV/c^2]
        self.steps_per_turn = 100
        self.time_per_turn = 1e-6 # [seconds]
        self.run_name = None # default is "PyOpal"
        self.exit_code = 0

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
        beam.beam_frequency = 1e-6/self.time_per_turn # MHz
        beam.number_of_slices = 10
        beam.number_of_particles = int(self.distribution_str.split()[0])
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
        drift.end_normal_x=0.0
        drift.end_normal_y=1.0
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
        self.ring.maximum_r = self.r0*100
        self.ring.is_closed = False

    def make_option(self):
        """Options enable setting of global control variables.

        No options are set by default. For a full list of variables see the
        Option docs.
        """
        self.option = pyopal.objects.option.Option()
        self.option.execute()

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
        track.max_steps = [self.max_steps]
        track.steps_per_turn = self.steps_per_turn
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
        self.track_run = run

    def make_element_iterable(self):
        """
        Return an iterable (e.g. list) of elements to append to the line

        By default, returns an empty list. User can overload this method.
        """
        return []

    def preprocess(self):
        """Perform any preprocessing steps just before the trackrun is executed

        This method can be overloaded with user required steps.
        """
        pass

    def postprocess(self):
        """Perform any postprocessing steps after the tracking is executed

        This method can be overloaded with user required steps.
        """
        pass

    def execute(self):
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
            self.preprocess()
            self.track_run.execute()
            self.postprocess()
        except:
            raise
        finally:
            print("Finished running in directory", os.getcwd())
            os.chdir(here)

    def execute_fork(self):
        """
        Set up and run a simulation in a fork of the current process. The 

        This method is memory safe - resources are only created in the
        forked process, which is destroyed when the process concludes. The
        downside is that resources (e.g. lattice objects, etc) are destroyed
        when the process concludes. If something is needed, overload the
        preprocess and postprocess routines to do anything just before or just
        after tracking.

        Returns the return code of the forked process, given by self.exit_code
        from the forked MinimalRunner. This can be used as a simple flag for
        comms from the child to the parent process (e.g. for testing purposes).
        For example, postprocess(self) can be overloaded to set an exit code.

        Tested in linux, I don't know about OSX. Unlikely to work in any
        Windows environment.
        """
        a_pid = os.fork()
        if a_pid == 0: # the child process
            self.execute()
            # hard exit returning exit_code - don't want to end up in any exit
            # handling stuff, just die ungracefully now the simulation has run
            os._exit(self.exit_code)
        else:
            retvalue = os.waitpid(a_pid, 0)[1]
        return retvalue

    distribution_str = """1
0.0 0.0 0.0 0.0 0.0 0.0
"""

    ring_error = "Failed to append ring to the line. Try running make_ring()"+\
    " before calling make_line()."

def main():
    runner = MinimalRunner()
    runner.execute_fork()

if __name__ == "__main__":
    main()
