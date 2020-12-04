#import PyOpal.parser
import PyOpal.ndgrid

print("Default constructor")
PyOpal.ndgrid.NDGrid()
print("Variable spacing")
grid = PyOpal.ndgrid.NDGrid.initialise_variable_spacing([[10, 11, 12], [2, 3, 4]])
print(grid.size(1))
print(grid.coord_vector(1))
print(grid.coord_vector(2))

print("Constant spacing")
PyOpal.ndgrid.NDGrid.initialise_fixed_spacing([10, 11, 12], [2, 3, 4], [0, 1, 2])
print(grid.size(2))
