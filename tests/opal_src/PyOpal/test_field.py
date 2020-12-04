import sys
import os
import matplotlib
import unittest

import pyopal
import pyopal.parser
import pyopal.field

class TestField(unittest.TestCase):
	def test_throws_on_no_lattice(self):
		try:
			field = pyopal.field.get_field_value(0.5, 5, 0, 0)
			self.assertTrue(false, "Should have thrown an exception")
		except RuntimeError:
			pass # Lattice not initialised - correctly threw error

	def test_initialise_from_opal_file(self):
		test_dir = os.path.split(__file__)[0]
		lattice_dir = os.path.join(test_dir, "dummy_lattice")
		os.chdir(lattice_dir)
		pyopal.parser.initialise_from_opal_file("dummy_lattice.opal")

		field = pyopal.field.get_field_value(0.5, 5, 0, 0)
		print("field      ", field)

		# Really I am testing the boost framework here; normally this is not
		# necessary, I am just doing it as a sort of experiment.
		try:
			field = pyopal.field.get_field_value("fish", 5, 0, 0)
			self.assertTrue(false, "Should have thrown an exception")
		except TypeError:
			pass # bad value - correctly threw an error

if __name__ == "__main__":
	unittest.main() # for running directly from the command line
