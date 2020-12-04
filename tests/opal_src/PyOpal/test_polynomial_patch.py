import PyOpal.ndgrid
import PyOpal.polynomial_patch


grid = PyOpal.ndgrid.NDGrid.initialise_fixed_spacing([2, 2], [1, 1], [0, 0])

values = [
    [0, 0], [0, 1], [1, 0], [1, 1]
]
patch = PyOpal.polynomial_patch.PolynomialPatch.initialise_from_solve_factory(grid, values, 1, 1)
print(patch)
print(patch.function([0,0]))
print(patch.function([0,1]))
print(patch.function([1,1]))
print(patch.function([0.5,0.5]))
