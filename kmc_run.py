import numpy as np
import matplotlib.pyplot as plt
import time,sys

sys.path.insert(0,'lib/')
from KMCEngine import KMCEngine
from lattice_models import lattice_model
from plotter import run_and_plot

lat = lattice_model("graphene")

cycles = 400
steps_per_cycle = 2000

ncells = [100,int(100*np.tan(np.pi/8.))]

jump, bond, penalty = 1.0, 3.0, 0.0
theta = np.pi / 8.0
force = np.r_[[np.cos(theta),np.sin(theta)]] * 1.5

line_slopes = [0.,np.tan(np.pi/8.),np.tan(np.pi/4.)]


""" Initialize object """
sim = KMCEngine(cell=lat.unit_cell(), basis=lat.unit_basis(), radius=lat.rad, \
                jump=jump, bond=bond, force=force, penalty=penalty)
sim.build(ncells[0],ncells[1])
print "Built %d sites" % sim.nsites

""" Make initial loop (overwritten if occupation fiel exists)"""
pos = sim.get_positions()
ncells = [100,int(100*np.tan(np.pi/8.))]
init_pos = np.r_[[sim.super_cell[0]-4.*lat.rad,sim.super_cell[1]/2.-4.*lat.rad]]
center =  pos[np.linalg.norm(pos-init_pos,axis=1).argmin()]
sel = np.linalg.norm(pos-center,axis=1) > 4.0 * lat.rad
occupation = np.zeros(sim.nsites,bool)
occupation[sel] = True
sim.set_occupations(occupation)

""" define some guide lines """
lines = []
for sl in line_slopes:
    lines.append([sl,init_pos])

""" Set up and run the simulation """
run_and_plot(sim,cycles=200, steps_per_cycle=5000,\
                    snapshots=2,lines=lines,\
                    image_file = None,\
                    occupation_file='data/occupation',\
                    cluster_file='data/cluster',\
                    traj_file='data/traj')
                    #image_file='figs/graph_zag_zag.pdf',\
