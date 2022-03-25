from triangular_lattice import *
from make_mesh import *

t = TriangularLattice(21)
t.build_periodic_lattice()
m = MakeMesh(t)
m.make_initial_configuration()
m.json_out("hexagonal.json", params={'kappa': 1.0, 'gamma': 0.25}, stress_free=True)
