import numpy as np 
from make_polygon import *

# Output file name
out = 'poly_test.json'

# We are making a circle or radius R = 10
# We need to make the perimeter of the circle. It is important that perimeter points are close enough to each other, i.e.
# to be within the cutoff distance
R = 10 
rcut = 1.0
Nperim = int(round(2*np.pi*R/rcut))
phi = np.linspace(0,2*np.pi,Nperim)

# Make the perimeter
x = R*np.cos(phi)
y = R*np.sin(phi)
r = np.stack((x,y)).T 

# Prepare region
p = Polyogonal(r)

# Set the number of cells in the region
N = 100

# Generate initial seed points
p.generate_seeds(N, rcut = rcut)

# Optimise their positions
p.optimise(eps = 1e-5, maxiter = 5000)

# Generate actual cells using Voronooi diagrams
p.generate_cells()

# Calculate cell properties (area, perimeter, etc.)
p.compute_cell_properties()

# Print JSON file
# Please chance p0 accordingly
p.write_json(out, max_area_scale=2.0, p0 = None)  # maxA0 is set to max_area_scale*meanA, where meanA is the mean area of all internal cells (i.e., outer cell excluded)
