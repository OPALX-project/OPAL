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

"""Test the Field Solver"""

import os
import unittest
import pyopal.objects.field_solver

class TestFieldSolver(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Test the Field Solver"""
    def setUp(self):
        """Set up some data"""
        self.fs = pyopal.objects.field_solver.FieldSolver()
        self.fs.mesh_size_x = 2
        self.fs.mesh_size_y = 3
        self.fs.mesh_size_z = 4

        self.fs.parallelise_fft_x = True
        self.fs.parallelise_fft_y = True
        self.fs.parallelise_fft_t = True

        self.fs.fft_boundary_x = "Fred"
        self.fs.fft_boundary_y = "Jim"
        self.fs.fft_boundary_z = "Bob"

        self.fs.greens_function = "Sara"

        self.bounding_box_increase = 5

        self.fs.geometry = "Lucy"
        self.fs.iterative_solver = "Sally"
        self.fs.interpolation = "Cara"

        self.fs.tolerance = 6
        self.fs.max_iterations = 7
        self.fs.preconditioner_mode = "Billy"
        self.fs.cutoff_radius = 8
        self.fs.alpha = 9
        self.fs.epsilon = 10

    def encapsulated_test_init(self):
        """Check I didn't make any typos"""
        my_fs = pyopal.objects.field_solver.FieldSolver()
        self.assertEqual(self.fs.mesh_size_x, 2)
        self.assertEqual(self.fs.mesh_size_y, 3)
        self.assertEqual(self.fs.mesh_size_z, 4)

        self.assertEqual(self.fs.parallelise_fft_x, True)
        self.assertEqual(self.fs.parallelise_fft_y, True)
        self.assertEqual(self.fs.parallelise_fft_t, True)

        self.assertEqual(self.fs.fft_boundary_x, "FRED")
        self.assertEqual(self.fs.fft_boundary_y, "JIM")
        self.assertEqual(self.fs.fft_boundary_z, "BOB")

        self.assertEqual(self.fs.greens_function, "SARA")

        self.assertEqual(self.bounding_box_increase, 5)

        self.assertEqual(self.fs.geometry, "LUCY")
        self.assertEqual(self.fs.iterative_solver, "SALLY")
        self.assertEqual(self.fs.interpolation, "CARA")

        self.assertEqual(self.fs.tolerance, 6)
        self.assertEqual(self.fs.max_iterations, 7)
        self.assertEqual(self.fs.preconditioner_mode, "BILLY")
        self.assertEqual(self.fs.cutoff_radius, 8)
        self.assertEqual(self.fs.alpha, 9)
        self.assertEqual(self.fs.epsilon, 10)

    def encapsulated_test_register(self):
        """Check we can register the beam"""
        my_beam = pyopal.objects.field_solver.FieldSolver()
        my_beam.type = "None"
        my_beam.register()


if __name__ == "__main__":
    unittest.main()
