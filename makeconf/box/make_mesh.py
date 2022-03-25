from cell import Cell
from vertex import Vertex
from random_box import RandomBox

from scipy.spatial import Voronoi
import numpy as np

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

import json
import copy


class MakeMesh:

    def __init__(self, Lx, Ly=None):
        self.Lx = Lx
        self.Ly = Ly if Ly != None else Lx

    def construct_random(self, points):
        self.random_box = RandomBox(len(points), self.Lx, self.Ly)
        self.random_box.set_points(points)
        self.random_box.reflect(3.0)
        all_points = np.vstack(
            (self.random_box.get_coords(), self.random_box.added_points))
        vor = Voronoi(all_points)
        self.mesh_vertices = []
        self.mesh_faces = []
        vert_map = {}
        vidx = 0
        for i in range(len(self.random_box.get_coords())):
            tile = vor.regions[vor.point_region[i]]
            tile_verts = []
            for v in tile:
                if not v in vert_map:
                    vert_map[v] = vidx
                    self.mesh_vertices.append(Vertex(vidx, vor.vertices[v]))
                    vidx += 1
                tile_verts.append(vert_map[v])
            c = Cell(i, vor.points[i])
            c.verts = copy.copy(tile_verts)
            self.mesh_faces.append(c)
        for c in self.mesh_faces:
            self.order_cell_vertices(c)
        # Build outer face
        c = Cell(i+1, np.array([0, 0]))
        c.outer = True
        xmin, xmax = -0.5*self.Lx, 0.5*self.Lx
        ymin, ymax = -0.5*self.Ly, 0.5*self.Ly
        for vi in range(len(self.mesh_vertices)):
            x, y = self.mesh_vertices[vi].r
            if abs(x - xmin) <= 1e-5 or abs(x - xmax) <= 1e-5 or abs(y - ymin) <= 1e-5 or abs(y - ymax) <= 1e-5:
                c.verts.append(vi)
                self.mesh_vertices[vi].boundary = True
        self.order_cell_vertices(c, True)
        self.mesh_faces.append(c)

    def cell_area(self, c):
        r = []
        for v in c.verts:
            r.append(self.mesh_vertices[v].r)
        r = np.array(r)
        i = np.arange(len(c.verts))
        j = np.roll(i, -1)
        return 0.5 * np.sum(r[i, 0] * r[j, 1] - r[j, 0] * r[i, 1])

    def cell_perim(self, c):
        r = []
        for v in c.verts:
            r.append(self.mesh_vertices[v].r)
        r = np.array(r)
        return np.sum(np.sqrt((r[1:, 0]-r[:-1, 0])**2 + (r[1:, 1]-r[:-1, 1])**2)) + np.sqrt((r[-1, 0]-r[0, 0])**2 + (r[-1, 1]-r[0, 1])**2)

    def order_cell_vertices(self, c, clockwise=False):
        rc = c.rc
        angles = []
        for v in c.verts:
            dx, dy = self.mesh_vertices[v].r - rc
            angles.append(np.arctan2(dy, dx))
        sidx = np.argsort(angles)
        if clockwise:
            sidx = list(reversed(sidx))
        verts = []
        for idx in sidx:
            verts.append(c.verts[idx])
        c.verts = copy.copy(verts)

    def set_cell_types(self, region, celltype):
        poly = Polygon(region)
        for c in self.mesh_faces:
            p = Point([c.rc[0], c.rc[1]])
            if poly.contains(p):
                c.type = celltype

    def set_vertex_types(self, region, verttype):
        poly = Polygon(region)
        for v in self.mesh_vertices:
            p = Point([v.r[0], v.r[1]])
            if poly.contains(p):
                v.type = verttype

    def mark_sides(self, sides):
        xmin, xmax = -0.5*self.Lx, 0.5*self.Lx
        ymin, ymax = -0.5*self.Ly, 0.5*self.Ly
        for vi in range(len(self.mesh_vertices)):
            x, y = self.mesh_vertices[vi].r
            if abs(x - xmin) <= 1e-5:
                if 'left' in sides:
                    self.mesh_vertices[vi].type = sides['left'][0]
                    if sides['left'][1]:
                        self.mesh_vertices[vi].constraint = 'x'
            elif abs(y - ymin) <= 1e-5:
                if 'bottom' in sides:
                    self.mesh_vertices[vi].type = sides['bottom'][0]
                    if sides['bottom'][1]:
                        self.mesh_vertices[vi].constraint = 'y'
            elif abs(x - xmax) <= 1e-5:
                if 'right' in sides:
                    self.mesh_vertices[vi].type = sides['right'][0]
                    if sides['right'][1]:
                        self.mesh_vertices[vi].constraint = 'x'
            elif abs(y - ymax) <= 1e-5:
                if 'top' in sides:
                    self.mesh_vertices[vi].type = sides['top'][0]
                    if sides['top'][1]:
                        self.mesh_vertices[vi].constraint = 'y'

    def set_elastic_parameters(self, params, stress_free=False):
        if stress_free:
            if 'lambda' in params:
                print(
                    'Warning! Using stress_free flag set. Lambda will be overwritten. ')
            if 'l0' in params:
                print('Warning! Using stress_free flag set. l0s will be overwritten. ')
        for c in self.mesh_faces:
            c.area = self.cell_area(c)
            c.perim = self.cell_perim(c)
            if 'kappa' in params:
                c.kappa = params['kappa']
            if 'gamma' in params:
                c.gamma = params['gamma']
            if 'lambda' in params:
                c.lam = params['lambda']
            if 'k' in params:
                c.k = params['k']
            if 'l0' in params and not stress_free:
                c.l0.append(l0)
            if stress_free:
                if c.gamma == None:
                    raise Exception(
                        'Gamma needs to be set before we can find lambda for the stress free configuration.')
                c.lam = c.gamma*c.perim
                N = len(c.verts)
                for i in range(N):
                    vi = self.mesh_vertices[c.verts[i]]
                    vj = self.mesh_vertices[c.verts[(i+1) % N]]
                    c.l0.append(
                        np.sqrt((vi.r[0] - vj.r[0])**2 + (vi.r[1] - vj.r[1])**2))

    def set_cell_myo(self, cell_myo=None, junction_myo=None):
        for c in self.mesh_faces:
            c.cell_myo = cell_myo if cell_myo != None else np.random.uniform(
                4, 8)
            if junction_myo != None:
                for i in range(len(c.verts)):
                    c.myo.append(np.random.uniform(0, 1))

    def json_out(self, fname):
        jsonData = {}
        jsonData["mesh"] = {}
        jsonData["mesh"]["vertices"] = []
        jsonData["mesh"]["l0"] = 1.0
        jsonData["mesh"]["box"] = {
            "periodic": False, "lx": self.Lx, "ly": self.Ly}
        for i in range(len(self.mesh_vertices)):
            vd = {}
            vd["id"] = i
            vd["r"] = self.mesh_vertices[i].r.tolist()
            vd["type"] = self.mesh_vertices[i].type
            vd["erased"] = False
            vd["boundary"] = self.mesh_vertices[i].boundary
            vd["constraint"] = self.mesh_vertices[i].constraint
            jsonData["mesh"]["vertices"].append(vd)
        jsonData["mesh"]["faces"] = []
        for c in self.mesh_faces:
            cd = {}
            cd["id"] = c.id
            cd["outer"] = c.outer
            cd["nsides"] = len(c.verts)
            cd["type"] = c.type
            cd["vertices"] = copy.copy(c.verts)
            cd["A0"] = self.cell_area(c) if c.area == None else c.area
            cd["P0"] = self.cell_perim(c) if c.perim == None else c.perim
            if c.kappa != None:
                cd['kappa'] = c.kappa
            if c.gamma != None:
                cd['gamma'] = c.gamma
            if c.lam != None:
                cd['lambda'] = c.lam
            if c.k != None:
                cd['k'] = c.k
            if len(c.l0) > 0:
                cd['l0'] = copy.copy(c.l0)
            if c.cell_myo != None:
                cd['cell_myo'] = c.cell_myo
            if len(c.myo) > 0:
                cd['myo'] = copy.copy(c.myo)
            jsonData["mesh"]["faces"].append(cd)
        with open(fname, 'w') as out:
            json.dump(jsonData, out, sort_keys=True, indent=4)
