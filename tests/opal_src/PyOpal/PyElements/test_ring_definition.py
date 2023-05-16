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

"""Test RingDefinition python implementation"""
import math
import os
import unittest
import pyopal.objects.encapsulated_test_case
import pyopal.objects.minimal_runner
import pyopal.objects.field
import pyopal.elements.local_cartesian_offset
import pyopal.elements.ring_definition

class RingDefinitionRunner(pyopal.objects.minimal_runner.MinimalRunner):
    def __init__(self):
        super().__init__()
        self.run_name = "ring_test"
        self.momentum = 0.1 # [GeV/c]
        self.mass = 0.93827208816 # [GeV/c^2]
        self.verbose = 0

    def make_ring(self):
        self.ring = pyopal.elements.ring_definition.RingDefinition()
        self.ring.set_opal_name("a_ring")
        self.ring.set_attributes(
            lattice_initial_r = 5.0,
            lattice_initial_phi = 45.0,
            lattice_initial_theta = 45.0,
            beam_initial_r = 5.5,
            beam_initial_phi = 45.0,
            beam_initial_theta = 45.0
        )

class RingDefinitionTest(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Test RingDefinition python implementation"""
    def setUp(self):
        """setup the ring"""
        self.ring_definition = pyopal.elements.ring_definition.RingDefinition()
        self.ring_definition.lattice_initial_r = 1.0
        self.ring_definition.lattice_initial_phi = 2.0
        self.ring_definition.lattice_initial_theta = 3.0
        self.ring_definition.beam_initial_r = 4.0
        self.ring_definition.beam_initial_phi = 5.0
        self.ring_definition.beam_initial_theta = 5.5
        self.ring_definition.beam_initial_pr = 6.0
        self.ring_definition.harmonic_number = 7.0
        self.ring_definition.symmetry = 8
        self.ring_definition.scale = 9.0
        self.ring_definition.rf_frequency = 10.0
        self.ring_definition.is_closed = True
        self.ring_definition.minimum_r = 12.0
        self.ring_definition.maximum_r = 13.0

    def test_init(self):
        """Check that we can initialise the ring"""
        self.assertEqual(self.ring_definition.lattice_initial_r, 1.0)
        self.assertEqual(self.ring_definition.lattice_initial_phi, 2.0)
        self.assertEqual(self.ring_definition.lattice_initial_theta, 3.0)
        self.assertEqual(self.ring_definition.beam_initial_r, 4.0)
        self.assertEqual(self.ring_definition.beam_initial_phi, 5.0)
        self.assertEqual(self.ring_definition.beam_initial_theta, 5.5)
        self.assertEqual(self.ring_definition.beam_initial_pr, 6.0)
        self.assertEqual(self.ring_definition.harmonic_number, 7.0)
        self.assertEqual(self.ring_definition.symmetry, 8)
        self.assertEqual(self.ring_definition.scale, 9.0)
        self.assertEqual(self.ring_definition.rf_frequency, 10.0)
        self.assertEqual(self.ring_definition.is_closed, True)
        self.assertEqual(self.ring_definition.minimum_r, 12.0)
        self.assertEqual(self.ring_definition.maximum_r, 13.0)

    def test_get_field_value(self):
        """
        Check get_field_value on the fields

        Field will always return outofbounds and 0 - we can only add elements
        using TrackRun, so we will test get_field_value then.
        """
        value = self.ring_definition.get_field_value(11.9, 0.0, 0.0, 0.0)
        self.assertTrue(value[0])
        for i in range(1, 7):
            self.assertEqual(value[i], 0.0)

    def encapsulated_test_integration(self):
        ring_def = RingDefinitionRunner()
        ring_def.execute()
        position = pyopal.objects.field.get_element_start_position(0)
        normal = pyopal.objects.field.get_element_start_normal(0)
        # radius is 5.0
        self.assertAlmostEqual(position[0]**2+position[1]**2, 25.0, 6)
        # start angle is 45 degrees
        self.assertAlmostEqual(math.atan2(position[1], position[0]), math.radians(45))
        # start angle is 90 degrees including lattice "theta" i.e. pointing towards
        # negative x
        self.assertAlmostEqual(normal[0], -1.0)
        self.assertAlmostEqual(normal[1], 0.0)

        track_orbit = os.path.join(ring_def.tmp_dir, "ring_test-trackOrbit.dat")
        with open(track_orbit) as fin:
            fin.readline()
            fin.readline()
            words = fin.readline().split()
        self.assertAlmostEqual(float(words[1]), 5.5/2**0.5) # x
        self.assertAlmostEqual(float(words[3]), 5.5/2**0.5) # y
        p = float(words[2]), float(words[4])
        self.assertAlmostEqual(p[0]+(p[0]**2+p[1]**2)**0.5, 0.0) # px
        self.assertAlmostEqual(float(p[1]), 0.0) #py



if __name__ == "__main__":
    unittest.main()
