#
#  Cell class
#
#  Cell list class
#
#

from copy import deepcopy

# Single cell of the CellList object
# it has its own index (label), position r (==(rx,ry)) and extent L (==(Lx,Ly))


class CLCell:

    def __init__(self, idx, r, L):
        self.idx = idx
        self.r = r
        self.L = L
        self.indices = []
        # Labels (indices) of the neighbouring cells
        self.neighbors = []

    def add_particle(self, idx):
        self.indices.append(idx)

    def printMe(self):
        print("I am cell " + str(self.idx))
        print("My position is: " + str(self.r))
        print("My particles are: " + str(self.indices))
        print("My neighbour cells are: " + str(self.neighbors))


class CellList2D:
    # Create my boxes
    def __init__(self, Lx, Ly=None, r_cut=1.0):
        self.Lx = Lx
        self.Ly = Ly if Ly != None else Lx
        self.r_cut = r_cut
        self.cell_indices = {}
        self.nx = int(self.Lx/r_cut)
        self.ny = int(self.Ly/r_cut)
        self.dx = self.Lx/float(self.nx)
        self.dy = self.Ly/float(self.ny)
        # total number of cells
        n_cell = self.nx*self.ny
        print("Created CellList with " + str(n_cell) + " squares.")
        # Cell list is a python list
        self.cell_list = [None for i in range(n_cell)]
        for i in range(self.nx):
            x = -0.5*self.Lx + float(i)*self.dx
            for j in range(self.ny):
                y = -0.5*self.Ly + float(j)*self.dy
                # Cell labeling scheme: for each x, do all y, and for all y, do all z
                idx = self.ny*i + j
                # Create new cell with index, position and size
                self.cell_list[idx] = CLCell(
                    idx, (x, y, 0.0), (self.dx, self.dy, 0.0))
                for ix in [-1, 0, 1]:
                    for iy in [-1, 0, 1]:
                        iix, iiy = i + ix, j + iy
                        if iix == self.nx:
                            iix = 0
                        elif iix < 0:
                            iix = self.nx - 1
                        if iiy == self.ny:
                            iiy = 0
                        elif iiy < 0:
                            iiy = self.ny - 1
                        self.cell_list[idx].neighbors.append(self.ny*iix + iiy)

    # Get the cell label for a given position vector v
    def get_cell_idx(self, rval):
        xmin, ymin = -0.5*self.Lx, -0.5*self.Ly
        i, j = int((rval[0]-xmin)/self.dx), int((rval[1]-ymin)/self.dy)
        cell_idx = self.ny*i + j
        return cell_idx

    # Add a particle to a cell: This means compute its cell index (if not given already)
    # Then add it with add_vertex of cell
    # Give the index to the list of cell indices: this particle is in this cell
    def add_particle(self, rval, idx, cell_idx=None):
        if cell_idx == None:
            cell_idx = self.get_cell_idx(rval)
        else:
            cell_idx = cell_index
        self.cell_list[cell_idx].add_particle(idx)
        self.cell_indices[idx] = cell_idx

    # Nuke the contents of the whole cell list

    def wipe(self):
        for cell in self.cell_list:
            cell.vertices = []

    # Get the neighbours of particle at position rval
    # This means looking at the neighbour boxes (including oneself), and copying their list of neighbours into one
    # single neighbour list
    def get_neighbours(self, rval):
        cell_index = self.get_cell_idx(rval)
        neighbors = []
        for idx in self.cell_list[cell_index].neighbors:
            neighbors.extend(deepcopy(self.cell_list[idx].indices))
        return neighbors

    def printMe(self):
        for cell in self.cell_list:
            cell.printMe()
