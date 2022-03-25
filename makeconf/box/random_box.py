import numpy as np
from vertex import Vertex
from CellList2D import CellList2D

# Randomly place paticles in a rectangular box


class RandomBox:

    def __init__(self, N, Lx, Ly=None, min_dist=1.0):
        self.N = N
        self.Lx = Lx
        self.Ly = Ly if Ly != None else Lx
        self.min_dist = 1.0
        self.cl = CellList2D(Lx, Ly, 2.5*min_dist)
        self.built = False
        self.max_attempts = 5*N

    def build_lattice(self):
        i = 0
        points = []
        attempt = 0
        while i < self.N:
            x, y = np.random.uniform(-0.5*(self.Lx-self.min_dist), 0.5*(self.Lx-self.min_dist)
                                     ), np.random.uniform(-0.5*(self.Ly-self.min_dist), 0.5*(self.Ly-self.min_dist))
            r = np.array([x, y])
            neighs = self.cl.get_neighbours(r)
            can_add = True
            for n in neighs:
                vn = points[n]
                rn = vn.r
                dr = rn - r
                if np.sqrt(np.dot(dr, dr)) < self.min_dist:
                    can_add = False
                    break
            if can_add:
                points.append(Vertex(i, r))
                self.cl.add_particle(r, i)
                i += 1
            if attempt > self.max_attempts:
                raise Exception('Failed to place points after '+str(attempt) +
                                ' attempts. Please reduce the cutoff distance.')
            else:
                attempt += 1
        self.built = True
        self.points = np.array(points)

    def set_points(self, pts):
        points = []
        i = 0
        for p in pts:
            points.append(Vertex(i, p))
            i += 1
        self.N = len(points)
        self.points = np.array(points)

    def get_coords(self):
        coords = []
        for p in self.points:
            coords.append(p.r)
        return np.array(coords)

    def __reflect_point(self, r, l):
        p, q = r
        a, b, c = l
        d = a*a + b*b
        return np.array([(p*(a*a - b*b) - 2*b*(a*q + c))/d, (q*(b*b-a*a)-2*a*(b*p+c))/d])

    def reflect(self, pad_dist=3.0):
        xmin, xmax = -0.5*self.Lx, 0.5*self.Lx
        ymin, ymax = -0.5*self.Ly, 0.5*self.Ly
        cutoff_dist = pad_dist*self.min_dist
        added_points = []
        for p in self.points:
            r = p.r
            # left
            if r[0] - xmin <= cutoff_dist:
                added_points.append(self.__reflect_point(r, [0, 1, xmax]))
            # right
            if xmax - r[0] <= cutoff_dist:
                added_points.append(self.__reflect_point(r, [0, 1, xmin]))
            # bottom
            if r[1] - ymin <= cutoff_dist:
                added_points.append(self.__reflect_point(r, [1, 0, ymax]))
            # top
            if ymax - r[1] <= cutoff_dist:
                added_points.append(self.__reflect_point(r, [1, 0, ymin]))
        self.added_points = np.array(added_points)
