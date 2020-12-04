import sys
import os
import matplotlib

import PyOpal
import PyOpal.parser
import PyOpal.field
import PyOpal.boost_field
import PyOpal.bind_field

print(PyOpal.__file__)

os.chdir("dummy_lattice")

try:
	field = PyOpal.boost_field.get_field_value(0.5, 5, 0, 0)
except RuntimeError:
	print("Lattice not initialised - correctly threw error")
try:
	field = PyOpal.field.get_field_value(0.5, 5, 0, 0)
except ValueError:
	print("Lattice not initialised - correctly threw error")
try:
	field = PyOpal.bind_field.get_field_value(0.5, 5, 0, 0)
except RuntimeError:
	print("Lattice not initialised - correctly threw error")

input("Press <CR> to initialise lattice")

PyOpal.parser.initialise_from_opal_file("dummy_lattice.opal")
input("Press <CR> to get some fields")

field = PyOpal.boost_field.get_field_value(0.5, 5, 0, 0)
print("boost_field", field)
field = PyOpal.field.get_field_value(0.5, 5, 0, 0)
print("field      ", field)
field = PyOpal.bind_field.get_field_value(0.5, 5, 0, 0)
print("bind_field ", field)

input("Press <CR> to check error handling (not a double)")
try:
	field = PyOpal.boost_field.get_field_value("fish", 5, 0, 0)
except TypeError:
	pass
try:
	field = PyOpal.bind_field.get_field_value("fish", 5, 0, 0)
except TypeError:
	pass
try:
	field = PyOpal.field.get_field_value("fish", 5, 0, 0)
except TypeError:
	pass

print("Finished okay")
