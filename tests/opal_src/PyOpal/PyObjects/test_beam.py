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

"""Check that beam can be set up okay"""

import unittest

import pyopal.objects.minimal_runner
import pyopal.objects.encapsulated_test_case
import pyopal.objects.beam

class BeamRunner(pyopal.objects.minimal_runner.MinimalRunner):
    def __init__(self):
        super().__init__()
        self.run_name = "ring_test"
        self.momentum = 0.1 # [GeV/c]
        self.mass = 0.93827208816 # [GeV/c^2]
        self.momentum_tolerance = 1e-6
        self.verbose = 2

    def make_beam(self):
        beam = pyopal.objects.beam.Beam()
        beam.set_opal_name("SuperBeam")
        beam.mass = self.mass
        beam.momentum = self.momentum
        beam.charge = 1.0
        beam.beam_frequency = 1e-6/self.time_per_turn # MHz
        beam.number_of_slices = 10
        beam.number_of_particles = int(self.distribution_str.split()[0])
        beam.momentum_tolerance = self.momentum_tolerance
        beam.register()
        self.beam = beam

    def set_distribution_momentum(self, momentum):
        momentum = (momentum-self.momentum)*1e9
        self.distribution_str = f"""1
0.0 0.0 0.0 {momentum} 0.0 0.0
"""

class TestBeam(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Very light beam test class"""
    def encapsulated_test_init(self):
        my_beam = pyopal.objects.beam.Beam()
        my_beam.set_attributes(
            particle = "proton",
            mass = 1,
            charge = 2,
            energy = 3,
            momentum = 4,
            gamma = 5,
            beam_current = 6,
            beam_frequency = 7,
            number_of_particles = 8,
            momentum_tolerance = 0.001,
        )
        self.assertEqual(my_beam.particle, "PROTON")
        self.assertAlmostEqual(my_beam.mass, 1)
        self.assertAlmostEqual(my_beam.charge, 2)
        self.assertAlmostEqual(my_beam.energy, 3)
        self.assertAlmostEqual(my_beam.momentum, 4)
        self.assertAlmostEqual(my_beam.gamma, 5)
        self.assertAlmostEqual(my_beam.beam_current, 6)
        self.assertAlmostEqual(my_beam.beam_frequency, 7)
        self.assertAlmostEqual(my_beam.number_of_particles, 8)
        self.assertAlmostEqual(my_beam.momentum_tolerance, 0.001)

    def encapsulated_test_bad_particle_type(self):
        """See what happens if we set a bad particle type"""
        my_beam = pyopal.objects.beam.Beam()
        my_beam.particle = "bad"
        # would not be terrible if this threw exception

    def encapsulated_test_register(self):
        """Check we can register the beam okay"""
        my_beam = pyopal.objects.beam.Beam()
        my_beam.register()

    def encapsulated_test_momentum_tolerance(self):
        runner = BeamRunner()

        exit_code = runner.execute_fork()
        self.assertEqual(exit_code, 0) # all okay

        runner.set_distribution_momentum(0.10001)
        runner.momentum_tolerance = 1e-6
        exit_code = runner.execute_fork()
        self.assertNotEqual(exit_code, 0) # fails tolerance check

        runner.momentum_tolerance = 1e-3
        exit_code = runner.execute_fork()
        self.assertEqual(exit_code, 0) # back inside tolerance

        runner.set_distribution_momentum(0.101)
        exit_code = runner.execute_fork()
        self.assertNotEqual(exit_code, 0) # outside tolerance

        runner.momentum_tolerance = 0 # disable checks
        exit_code = runner.execute_fork()
        self.assertEqual(exit_code, 0) # no check, all okay


if __name__ == "__main__":
    unittest.main()
