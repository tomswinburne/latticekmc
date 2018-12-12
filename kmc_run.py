import numpy as np
import matplotlib.pyplot as plt
import time,sys

sys.path.insert(0,'lib/')
from KMCEngine import KMCEngine
from lattice_models import lattice_model
from simulations import simulator

cycles = 400
steps_per_cycle = 2000

jump = 1.0
bond = 3.0 # energy per bond
penalty = 0.

theta = 0.*np.pi / 8.0
force = 1.5 * np.r_[[np.cos(theta),np.sin(theta)]]

line_slopes = [0.,np.tan(-np.pi/4.),np.tan(np.pi/4.)]

lat = lattice_model("fcc")
ncells = [100,int(100*np.tan(np.pi/8.))]

""" Initialize object """
sim = KMCEngine(cell=lat.unit_cell(), basis=lat.unit_basis(), radius=lat.rad,\
                jump=jump, bond=bond, force=force, penalty=penalty)
sim.build(ncells[0],ncells[1])

print "%d sites" % sim.nsites


""" Get positions """
pos = sim.get_positions()

""" Make initial loop (overwritten if occupation fiel exists)"""
init_pos = np.r_[[sim.super_cell[0]-4.*lat.rad,sim.super_cell[1]/2.-4.*lat.rad]]
center =  pos[np.linalg.norm(pos-init_pos,axis=1).argmin()]
sel = np.linalg.norm(pos-center,axis=1) > 4.0 * lat.rad
occupation = np.zeros(sim.nsites,bool)
occupation[sel] = True

""" set occupation """
sim.set_occupations(occupation)

""" define some guide lines """
lines = []
for sl in line_slopes:
    lines.append([sl,init_pos])


traj_simulation = simulator(sim,cycles=2000, steps_per_cycle=5000,\
                    snapshots=10, file=None, lines=lines,\
                    occupation_file='data/occupation',\
                    cluster_file='data/cluster',traj_file='data/traj')
traj_simulation.run()


""" Run the simulation """
#zoetrope_simulation(sim,nframes=12,yf=3,steps=20000,file=None,lines=None)
#traj_simulation(sim, cycles=cycles, steps_per_cycle=steps_per_cycle,\
#                    file="hcp.pdf", lines=lines)

exit()
