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
        self.phase = self.make_time_dependence("PY_PHASE", 0.0, 0.0)
        self.voltage = self.make_time_dependence("PY_VOLTAGE", 2.0, 0.0)
        self.frequency = self.make_time_dependence("PY_FREQUENCY", 3.0, 4.0)
        self.rf = pyopal.elements.variable_rf_cavity.VariableRFCavity()
        self.rf.phase_model = "PY_PHASE"
        self.rf.amplitude_model = "PY_VOLTAGE"
        self.rf.frequency_model = "PY_FREQUENCY"
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
        """Check that field value is finite and time-dependent."""
        values = []
        for it in range(-500, 1501, 100):
            t = it*1e-9
            e_z = self.rf.get_field_value(0.0, 0.0, 0.0, t)[6]
            self.assertTrue(math.isfinite(e_z))
            values.append(e_z)
        self.assertGreater(max(values)-min(values), 1e-6)


if __name__ == "__main__":
    unittest.main()
