from random_lattice import *
from make_mesh import *
import os 

seed = int(os.urandom(4).hex(), base=16)  # generate some large number as seed, it will we stored in the JSON file

t = RandomLattice(240, 30, 27.712812921102035, 1.0, seed = seed)   # change second and third arguments to make a larger system
t.build_periodic_lattice()
m = MakeMesh(t)
m.make_initial_configuration()
#m.json_out('random_test.json')
#m.json_out('random_test.json', {'kappa': 1.0, 'gamma': 0.25}, stress_free = True)
m.json_out('random_test.json', {'A0': 3.4641016151377544, 'kappa': 1.0, 'gamma': 0.25, 'lambda': 1.736})
