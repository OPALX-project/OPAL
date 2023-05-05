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
Test variable RF cavity in python
"""

import math
import unittest
import pyopal.objects.minimal_runner

import pyopal.objects.option
import pyopal.elements.polynomial_time_dependence
import pyopal.elements.variable_rf_cavity


class TestVariableRFCavity(unittest.TestCase):
    """
    Test variable RF cavity in python
    """
    @classmethod
    def make_time_dependence(cls, name, pol0, pol1):
        """Make a time dependence"""
        ptd = pyopal.elements.polynomial_time_dependence.PolynomialTimeDependence()
        ptd.p0 = pol0
        ptd.p1 = pol1
        ptd.set_opal_name(name)
        ptd.update()
        return ptd

    def setUp(self):
        """Set up the cavity"""
        self.phase = self.make_time_dependence("phase", 0.0, 0.0)
        self.voltage = self.make_time_dependence("voltage", 2.0, 0.0)
        self.frequency = self.make_time_dependence("frequency", 3.0, 4.0)
        self.rf = pyopal.elements.variable_rf_cavity.VariableRFCavity()
        self.rf.phase_model = "phase"
        self.rf.amplitude_model = "voltage"
        self.rf.frequency_model = "frequency"
        self.rf.height = 0.5
        self.rf.width = 0.7
        self.rf.length = 0.8

    def test_bad_phase_model(self):
        """Check that rf cavity setup fails if no TD set"""
        self.rf.phase_model = "dummy"
        with self.assertRaises(RuntimeError):
            self.rf.update_time_dependence()

    def test_bad_amp_model(self):
        """Check that rf cavity setup fails if no TD set"""
        self.rf.amplitude_model = "dummy"
        with self.assertRaises(RuntimeError):
            self.rf.update_time_dependence()

    def test_bad_freq_model(self):
        """Check that rf cavity setup fails if no TD set"""
        self.rf.frequency_model = "dummy"
        with self.assertRaises(RuntimeError):
            self.rf.update_time_dependence()

    def test_bb(self):
        """Check that height and width are set okay"""
        self.assertFalse(self.rf.get_field_value(0.0, 0.0, 0.79, 0.0)[0])
        self.assertTrue(self.rf.get_field_value(0.0, 0.0, 0.81, 0.0)[0])
        self.assertFalse(self.rf.get_field_value(0.0, 0.24, 0.4, 0.0)[0])
        self.assertTrue(self.rf.get_field_value(0.0, 0.26, 0.4, 0.0)[0])
        self.assertFalse(self.rf.get_field_value(0.34, 0.0, 0.4, 0.0)[0])
        self.assertTrue(self.rf.get_field_value(0.36, 0.0, 0.4, 0.0)[0])

    def test_field(self):
        """Check that field value returns okay"""
        for it in range(10):
            t = it*1e-3
            e_z = self.rf.get_field_value(0.0, 0.0, 0.0, t)[6]
            freq = self.frequency.function(t)
            v_0 = self.voltage.function(t)
            ez_test = v_0*math.sin(2*math.pi*t*1e-3*freq) # NOTE the 1e-3 factor
            self.assertAlmostEqual(e_z, ez_test)


if __name__ == "__main__":
    unittest.main()
