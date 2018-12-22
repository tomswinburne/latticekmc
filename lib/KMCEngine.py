#!/usr/bin/env python
import os,sys,time
from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

def max_neighbours(R,basis,cell):
    # Find maximum number of sites in circle of radius R
    # To scan unit cell, need to make a square of side > 2R+cell_side
    sites = basis.copy()
    reps = [int(np.ceil(R/cell[0])),int(np.ceil(R/cell[1]))]
    i=0
    for rx in np.arange(-reps[0],reps[0]+1):
        for ry in np.arange(-reps[1],reps[1]+1):
            if rx==0 and ry==0:
                continue
            temp = basis.copy()
            temp[:,0] += rx*cell[0]
            temp[:,1] += ry*cell[1]
            sites = np.vstack((sites,temp))
            i+=1
    # go through basis atoms and take largest nn count
    nnc = 0
    for b in basis:
        nnc = max((np.linalg.norm(sites-b,axis=1)<R).sum()-1,nnc)
    return nnc



class KMCEngine:
    def __init__(self,cell=np.ones(2), \
                      basis=np.zeros((1,2)), radius=1.1, \
                      barriers=np.ones((5,5)), energies=-np.arange(5),\
                      force = np.zeros(2),name="lib/c++/libkmcsim.so"):

        self.sim = c_void_p()
        self.kmclib = cdll.LoadLibrary(name)

        self.kmclib.open_sim.argtypes = [\
                ndpointer(c_double, flags="C_CONTIGUOUS"),\
                ndpointer(c_double, flags="C_CONTIGUOUS"), c_uint, c_double,\
                ndpointer(c_double, flags="C_CONTIGUOUS"),\
                c_uint, ndpointer(c_double, flags="C_CONTIGUOUS"), c_void_p]

        self.kmclib.open_sim.restype = None

        self.kmclib.set_occupations.argtypes = [c_void_p,\
                                    ndpointer(c_bool, flags="C_CONTIGUOUS")]
        self.kmclib.set_occupations.restypes = None

        self.kmclib.get_positions.argtypes = [c_void_p,\
                                    ndpointer(c_double, flags="C_CONTIGUOUS")]
        self.kmclib.get_positions.restypes = None

        self.kmclib.get_neigh_count.argtypes = [c_void_p,\
                                    ndpointer(c_uint32, flags="C_CONTIGUOUS")]
        self.kmclib.get_neigh_count.restypes = None

        self.kmclib.get_occupations.argtypes = [c_void_p,\
                                    ndpointer(c_bool, flags="C_CONTIGUOUS")]
        self.kmclib.get_occupations.restypes = None

        self.kmclib.build.argtypes = [c_void_p, c_int, c_int]
        self.kmclib.build.restype = None

        self.kmclib.nsites.argtypes = [c_void_p]
        self.kmclib.nsites.restype = c_uint32

        self.kmclib.get_time.argtypes = [c_void_p]
        self.kmclib.get_time.restype = c_double

        self.kmclib.run.argtypes = [c_void_p, c_uint32, c_bool]
        self.kmclib.run.restype = None

        ebar = np.hstack((barriers.flatten(),energies.flatten()))
        self.kmclib.open_sim(cell.flatten(), basis.flatten(), basis.shape[0], \
                        radius,ebar, energies.shape[0], force, byref(self.sim))

        self.built = False
        self.now = time.time()
        self.unit_cell = cell
        self.unit_basis = basis
        self.super_cell = np.zeros(2)

        self.force = force
        self.barriers = barriers
        self.energies = energies
        self.radius = radius

    def build(self,xsize=10,ysize=None):
        if ysize is None:
            ysize = xsize
        self.kmclib.build(self.sim,xsize,ysize)
        self.nsites = self.kmclib.nsites(self.sim)
        self.start_time = 0.
        self.time = 0.
        self.super_cell[0] = self.unit_cell[0] * int(xsize/self.unit_cell[0])
        self.super_cell[1] = self.unit_cell[1] * int(ysize/self.unit_cell[1])
        self.built = True

    def get_time(self):
        self.time = self.kmclib.get_time(self.sim)
        return self.time + self.start_time

    def get_positions(self):
        if not self.built:
            exit(-1)
        data = np.zeros(2*self.nsites)
        self.kmclib.get_positions(self.sim,data)
        return data.reshape((-1,2))

    def get_neigh_count(self):
        if not self.built:
            exit(-1)
        data = np.zeros(self.nsites,np.uint32)
        self.kmclib.get_neigh_count(self.sim,data)
        return data.reshape((-1,1))

    def get_occupations(self):
        if not self.built:
            exit(-1)
        data = np.zeros(self.nsites,bool)
        self.kmclib.get_occupations(self.sim,data)
        return data.reshape((-1,1))

    def set_occupations(self,occv):
        self.kmclib.set_occupations(self.sim,occv)

    def run(self,steps=1,restart=True,seed=100):
        self.kmclib.run(self.sim,steps,restart,seed)
        if restart:
            print "Execution time: %2.2fs" % (time.time()-self.now)
            self.now = time.time()
