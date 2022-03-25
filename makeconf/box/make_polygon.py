from CellList2D import CellList2D
from vertex import Vertex
from cell import Cell

from scipy.spatial import Voronoi
import numpy as np

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from operator import itemgetter

import json
import copy

class Polyogonal:

    def __init__(self, boundary_points):
        self.seeds = boundary_points
        self.poly = Polygon(boundary_points)
        self.Nboundary = self.seeds.shape[0]
        self.has_points = False
    
    def generate_seeds(self, N, rcut = 2.0, max_attempts = 100):
        minx, miny, maxx, maxy = self.poly.bounds
        self.lx = maxx - minx 
        self.ly = maxy - miny
        self.points = [] 
        self.cl = CellList2D(1.5*self.lx, 1.5*self.ly, 2.5*rcut)
        for i in range(len(self.seeds)):
            self.cl.add_particle(self.seeds[i], i)
        offset = len(self.seeds)
        i = 0
        attempt = 0
        while i < N:
            x = np.random.uniform(minx,maxx)
            y = np.random.uniform(miny,maxy)
            if self.poly.contains(Point([x,y])):
                neighs = self.cl.get_neighbours([x,y])
                canadd = True
                for n in neighs:
                    xi, yi = self.seeds[n]
                    dl = np.sqrt((xi-x)**2 + (yi-y)**2)
                    if dl < rcut:
                        canadd = False
                        break 
                if canadd:
                    self.seeds = np.append(self.seeds,[[x,y]], axis=0)
                    self.cl.add_particle([x,y], i+offset)
                    i += 1
                    attempt = 0
                else:
                    attempt += 1
                if attempt > max_attempts:
                    raise RuntimeError('Maximum attempts reached.')
        self.has_points = True
    
    def optimise(self, eps = 1e-5, maxiter = 1000):
        if not self.has_points:
            raise RuntimeError('Initial seed points need to be generated with the generate_seeds function.')
        oldpos = np.copy(self.seeds[self.Nboundary:,:])
        newpos = np.zeros_like(oldpos)
        iteration = 0
        while iteration < maxiter:
            vor = Voronoi(self.seeds)
            i = 0
            for pr in vor.point_region[self.Nboundary:]:
                poly = Polygon(vor.vertices[vor.regions[pr]])
                x, y = poly.centroid.xy
                newpos[i,:] = [x[0],y[0]]
                i += 1
            dl = np.sqrt((newpos[:,0] - oldpos[:,0])**2 + (newpos[:,1] - oldpos[:,1])**2)
            if np.all(dl <= eps):
                break 
            else:
                print('Iteration : {:d}   eps = {:.6e}.'.format(iteration, np.max(dl)))
                self.seeds[self.Nboundary:,:] = newpos
                oldpos = np.copy(newpos)
                iteration += 1
        else:
            print('Optimisation did not converge.')
            

    def generate_cells(self):
        if not self.has_points:
            raise RuntimeError('Seed points need to be generated with the generate_seeds function.')
        vor = Voronoi(self.seeds)
        cells = []
        vertex_map = {}
        vidx = 0
        vertices = []
        for pr in vor.point_region[self.Nboundary:]:
            for v in vor.regions[pr]:
                if not v in vertex_map:
                    vertices.append(v)
                    vertex_map[v] = vidx 
                    vidx += 1
        boundary_vertices = []
        for pr in vor.point_region[:self.Nboundary]:
            for v in vor.regions[pr]:
                if v in vertices and not v in boundary_vertices:
                    boundary_vertices.append(v)
        vertex_map_list = list(map(lambda x : x[0], sorted(vertex_map.items(), key=itemgetter(1))))
        mapped_vertices = vor.vertices[vertex_map_list]
        self.vorocells = []
        for pr in vor.point_region[self.Nboundary:]:
            cell = []
            for v in vor.regions[pr]:
                cell.append(vertex_map[v])
            poly = Polygon(mapped_vertices[cell,:])
            xc, yc = poly.centroid.xy
            theta = np.arctan2(mapped_vertices[cell,1]-yc[0], mapped_vertices[cell,0]-xc[0])
            cidx = np.argsort(theta)
            self.vorocells.append(np.array(cell)[cidx].tolist())
        theta = np.arctan2(vor.vertices[boundary_vertices,1],vor.vertices[boundary_vertices,0])
        bidx = np.argsort(theta)
        cell = []
        for bi in bidx[::-1]:
            cell.append(vertex_map[boundary_vertices[bi]])
        self.boundary_vertices = []
        for b in boundary_vertices:
            self.boundary_vertices.append(vertex_map[b])
        self.vorocells.append(cell)
        self.vertices = []
        for i in range(mapped_vertices.shape[0]):
            self.vertices.append(Vertex(i,mapped_vertices[i,:]))
            if i in self.boundary_vertices:
                self.vertices[-1].boundary = True

    def compute_cell_properties(self):
        self.cells = []
        self.mean_area = 0.0
        for i in range(len(self.vorocells)):
            cell = self.vorocells[i]
            verts = []
            for v in cell:
                verts.append(self.vertices[v].r)
            poly = Polygon(verts)
            xc, yc = poly.centroid.xy
            c = Cell(i, [xc[0], yc[0]])
            c.verts = cell 
            c.area = poly.area
            c.perim = poly.length
            self.cells.append(c) 
            if i < len(self.vorocells) - 1:  # exclude outer cell
                self.mean_area += c.area
        self.cells[-1].outer = True
        self.cells[-1].area = -self.cells[-1].area  # outer cell has negative area by definition
        self.mean_area /= len(self.vorocells) - 1


    def write_json(self, fname, max_area_scale = 1.5, p0 = None):
        jsonData = {}
        jsonData["mesh"] = {}
        jsonData["mesh"]["vertices"] = []
        jsonData["mesh"]["box"] = {"periodic": False, "lx": self.lx, "ly": self.ly}
        for v in self.vertices:
            vd = {}
            vd["id"] = v.id
            vd["r"] = v.r.tolist()
            vd["type"] = v.type
            vd["erased"] = False 
            vd["boundary"] = v.boundary
            vd["constraint"] = v.constraint
            jsonData["mesh"]["vertices"].append(vd)
        jsonData["mesh"]["faces"] = []
        for c in self.cells:
            cd = {}
            cd["id"] = c.id
            cd["outer"] = c.outer
            cd["nsides"] = len(c.verts)
            cd["type"] = c.type
            cd["vertices"] = c.verts
            cd["A0"] = c.area
            if p0 == None:
                cd["P0"] = c.perim
            else:
                P0 = p0*np.sqrt(cd["A0"])
                cd["P0"] = P0
            cd["maxA0"] = max_area_scale*self.mean_area
            jsonData["mesh"]["faces"].append(cd)
        with open(fname,'w') as out:
            json.dump(jsonData, out, sort_keys = True, indent = 4) 
            



        
