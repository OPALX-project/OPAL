
import os
import unittest
import pyopal.objects.minimal_runner

import pyopal.elements.probe


class ProbeRunner(pyopal.objects.minimal_runner.MinimalRunner):
    def __init__(self):
        """Set up the test"""
        super().__init__()
        self.probe = None

    def make_probe(self):
        self.probe = pyopal.elements.probe.Probe()
        self.probe.x_start = 0.0
        self.probe.x_end = self.r0*2*1000
        self.probe.y_start = 1e-3
        self.probe.y_end = 1e-3
        self.probe.output_filename = "test_probe"
        self.probe.set_opal_name("my_probe")
        print(self.probe.get_opal_name())

    def make_element_iterable(self):
        self.make_probe()
        return [self.probe]


    def make_track(self):
        """Make a track object.

        The track object handles the interface between tracking and the beam
        elements.
        """
        track = pyopal.objects.track.Track()
        track.line = "test_line"
        track.beam = "SuperBeam"
        track.steps_per_turn = 1000
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



class TestProbe(unittest.TestCase):
    def test_one(self):
        probe_runner = ProbeRunner()
        probe_runner.run_one()
        track_file = os.path.join(probe_runner.tmp_dir, "PyOpal-trackOrbit.dat")
        with open(track_file) as a_file:
            print(a_file.read())

if __name__ == "__main__":
    unittest.main()