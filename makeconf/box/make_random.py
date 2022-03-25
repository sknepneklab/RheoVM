from voro_box import VoroBox
from make_mesh import MakeMesh
from shapely.geometry.polygon import Polygon
import numpy as np

N = 1100
Lx = 50
Ly = 50
minimal_interparticle_distance = 1.0
max_iterations = 2000
tolerance = 5e-5
optimised_coords = 'seeds.dat'
from_file = False  # Set this to True to build from a relaxed configuration
output = 'test.json'

# Optimise random cell packing
if not from_file:
    v = VoroBox(N, Lx, Ly, minimal_interparticle_distance)
    v.optimise(tol=tolerance, max_iter=max_iterations)
    np.savetxt(optimised_coords, v.random_box.get_coords())

# Build the mesh based on the optimal packin
m = MakeMesh(Lx, Ly)
if from_file:
    seeds = np.loadtxt(optimised_coords)
else:
    seeds = v.random_box.get_coords()

m.construct_random(seeds)
#m.mark_sides({'left' : ('left', True), 'bottom': ('bottom', False), 'right': ('right', True), 'top': ('top', False)})
m.mark_sides({'left': ('left', True), 'right': ('right', True)})

m.set_elastic_parameters({'kappa': 2.0, 'gamma': 0.25,
                         'lambda': 1.0, 'k': 0.5, 'l0': 1.0}, stress_free=True)
m.set_cell_myo(cell_myo=6.0)

# Set active region
active_region = [[-15, -20], [15, -20], [15, 0], [-15, 0]]
m.set_cell_types(active_region, 'active')
m.json_out(output)
