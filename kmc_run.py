#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import time,sys

sys.path.insert(0,'lib/')
from KMCEngine import KMCEngine, max_neighbours
from lattice_models import lattice_model
from plotter import run_and_plot

"""
TODO:
jump / connectivity histogram to investigate this strange reversal behaviour
custom jump matrix (hypermatrix with nnn) and energy vector (matrix with nnn)
"""


""" define lattice sites """

# simple 2d fcc
# We require orthorhombic unit cell so use vector not matrix (TODO: biclinic)
cell = np.r_[[1.,1.]]

# origin and fcc site
nbasis = 2
basis = np.zeros((nbasis,2))
basis[0,:] = np.r_[[0.,0.]]
basis[1,:] = np.r_[[.5,.5]]

# Number of simulation cells
ncells = [100,int(100*np.tan(np.pi/8.))]



""" define energy landscape """

# cut off for jump vectors
radius = 1.1/np.sqrt(2.)

# find maximum number of neighbours possible
max_nn = max_neighbours(radius,basis,cell)

# Define state energy vector and jump barrier matrix
# here we give simple bond counting + constant jump examples
bondE = 4.0 # in kT
jumpE = 1.0 # in kT
strength = .75 # in kT / length
theta = 0. # np.pi / 8.0


stateEn = -bondE * np.linspace(0.,1.*max_nn,max_nn+1,endpoint=True)

stateEn[0] = 6.0
print stateEn
jumpMat = jumpE * np.ones((max_nn+1,max_nn+1))
# needs to be symmetric for detailed balance!
#jumpMat = (jumpMat + jumpMat.T) / 2.

# Electomigration force
force = strength * np.r_[[np.cos(theta),np.sin(theta)]]

""" Initialize KMC object """
sim = KMCEngine(cell, basis, radius, jumpMat, stateEn, force)
sim.build(ncells[0],ncells[1])
print "Built %d sites" % sim.nsites

""" Make initial loop (overwritten if occupation file exists)"""
loop_radius = 4.*radius
pos = sim.get_positions()
init_pos = sim.super_cell/2.-loop_radius*np.ones(2)
init_pos[0] += sim.super_cell[0]/2.
init_pos =  pos[np.linalg.norm(pos-init_pos,axis=1).argmin()]
sel = np.linalg.norm(pos-init_pos,axis=1) > loop_radius
occupation = np.zeros(sim.nsites,bool)
occupation[sel] = True
sim.set_occupations(occupation)

""" simulation parameters """
cycles = 1000
steps_per_cycle = 10000
snapshotsevery = 500


""" define some guide lines """
lines = []
line_slopes = [0.,theta,np.tan(2.*theta)]
for sl in line_slopes:
    lines.append([sl,init_pos])


""" Set up and run the simulation """
run_and_plot(sim,cycles=cycles, steps_per_cycle=steps_per_cycle,\
                    snapshotsevery=snapshotsevery,lines=lines,\
                    occupation_file='data/temp/occupation',\
                    cluster_file='data/temp/cluster',\
                    traj_file='data/temp/traj',\
                    image_file = None)
                    #image_file='figs/graph_zag_zag.pdf',\
