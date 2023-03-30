# Import sys and os for data output paths
import sys as s
import os
import numpy as np

# Import rheovm module for simulatins
# Plase make sure to have typed in your shell export PYTHONPATH=$PYTHONPATH:/path/to/the/rheovm/build ( you need to do this only once)
from rheovm import *

arglist = s.argv
exponent = int(arglist[1])
num_period = int(arglist[2])

# This is where we start a simulation

tissue = Tissue()           # initialise mesh 
sys = System(tissue)        # base object for the system (mesh + parameters)
f = Force(sys)              # handles all types of forces  (all forces acting on a vertex)
integ = Integrate(sys,f,0)  # handles all integrators
t = Topology(sys, f)        # handles all topology changes (T1, division, ingression)
d = Dump(sys, f)            # handles all data output 
sim = Simulation(sys, integ, f, t)  # simulation object

sys.read_input('hexagons_L15.json')           # read input configuration
f.set_update(True)   # Keep track of all force components



#----------- Forces 
f.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
f.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P
f.add('self-propulsion') #add self-propulsion

# All force terms will be using paramters defines based on cell type (default is that parateters are given in the JSON file)
f.set_flag('area', 'use_cell_type')
f.set_flag('perimeter', 'use_cell_type')

f.set_params('area', 'passive', {'kappa' : 1.0})
f.set_params('perimeter', 'passive',  {'gamma': 0.25, 'lambda': 0.8})
f.set_compute_stress(True)

dt = 0.01

#----------- Integrators 
integ.add('brownian')    # add Brownian integrator that handles all vercecx movements
t.set_params({'min_edge_len': 0.2, 'new_edge_len': 0.22, 'myosin': 0.5}) # /!\ need to set new edges myosin value
integ.set_dt(dt)       # all integrators have the same time step size
integ.set_params('brownian', {'T' : 0.0})
integ.set_params('brownian',{'gamma' :0.1})
d.set_sfc(1.0)


gamma = 0.001*0.01*0.01   #magnitude of the shear
#num_period = 20
T = 2**exponent/4                #length of the period
omega = 2*np.pi/T
pts_period = int(T/dt)   #number of steps per period

np.savetxt('num_period_t_{:04.8f}.dat'.format(T), np.array([num_period]), fmt='%d')

if os.path.exists('stress_t_{:04.8f}.dat'.format(T)):
    os.remove('stress_t_{:04.8f}.dat'.format(T))

for i in range(0,num_period*pts_period):

    shear = gamma*np.sin(omega*dt*i) 
    tissue.transform(1,shear,0,1,undo=True)
    sim.run(1)

    #dump the stress 25 times in every cycle
    if  (i%(int(pts_period/25)) == 0):
        d.dump_stress('stress_t_{:04.8f}.dat'.format(T),True)

