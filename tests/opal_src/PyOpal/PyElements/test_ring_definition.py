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
import os
import unittest
import pyopal.objects.minimal_runner
import pyopal.elements.local_cartesian_offset
import pyopal.elements.ring_definition

class RingDefinitionTest(unittest.TestCase):
    """Test RingDefinition python implementation"""
    def setUp(self):
        """setup the ring"""
        self.ring_definition = pyopal.elements.ring_definition.RingDefinition()
        self.ring_definition.lattice_initial_r = 1.0
        self.ring_definition.lattice_initial_phi = 2.0
        self.ring_definition.lattice_initial_theta = 3.0
        self.ring_definition.beam_initial_r = 4.0
        self.ring_definition.beam_initial_phi = 5.0
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

if __name__ == "__main__":
    unittest.main()
