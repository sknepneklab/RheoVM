from random_box import RandomBox
from scipy.spatial import Voronoi
import numpy as np


class VoroBox:

    def __init__(self, N, Lx, Ly=None, min_dist=1.0):
        self.Lx = Lx
        self.Ly = Ly if Ly != None else Lx
        self.random_box = RandomBox(N, Lx, Ly, min_dist)
        self.random_box.build_lattice()

    def compute_centre(self, r):
        xc, yc = np.mean(r[:, 0]), np.mean(r[:, 1])
        return np.array([xc, yc])

    def optimise(self, tol=1e-6, max_iter=100):
        old_seeds = self.random_box.get_coords()
        iteration = 0
        while iteration < max_iter:
            self.random_box.reflect(3.0)
            inner_points = []
            for p in self.random_box.points:
                inner_points.append(p.r)
            all_points = np.vstack(
                (inner_points, self.random_box.added_points))
            new_seeds = []
            vor = Voronoi(all_points)
            for i in range(len(inner_points)):
                tile = vor.regions[vor.point_region[i]]
                r = vor.vertices[tile]
                new_seeds.append(self.compute_centre(r))
            new_seeds = np.array(new_seeds)
            dr = new_seeds - old_seeds
            dr2 = np.sqrt(dr[:, 0]**2 + dr[:, 1]**2)
            self.random_box.set_points(new_seeds)
            if np.max(dr2) < tol:
                break
            else:
                print('Iteration : {: 5d}, maximum centre displacement : {:.6e}'.format(
                    iteration, np.max(dr2)))
                old_seeds = np.copy(new_seeds)
                iteration += 1
        else:
            print(
                'Build did no converge after {:5d} iterations.'.format(max_iter))
