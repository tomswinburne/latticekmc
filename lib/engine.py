#!/usr/bin/env python
import os,sys,time
from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

class KMCEngine:
    def __init__(self,name="c++/libsim.so",\
            cell=np.ones(2),basis=np.zeros((1,2)),radius=1.,\
            jump=1.0, bond=4.0,force=np.zeros(2),penalty=1.):

        self.sim = c_void_p()
        self.kmclib = cdll.LoadLibrary(name)
        self.kmclib.open_sim.argtypes = [\
                ndpointer(c_double, flags="C_CONTIGUOUS"),\
                ndpointer(c_double, flags="C_CONTIGUOUS"),\
                c_int, c_double, c_double, c_double,\
                ndpointer(c_double, flags="C_CONTIGUOUS"),\
                c_double, c_void_p]
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

        self.kmclib.open_sim(cell,basis.flatten(),basis.shape[0],radius,\
                    jump, bond, force, penalty, byref(self.sim))
        print jump
        self.built = False
        self.now = time.time()



    def build(self,xsize=10,ysize=None):
        if ysize is None:
            ysize = xsize
        self.kmclib.build(self.sim,xsize,ysize)
        self.nsites = self.kmclib.nsites(self.sim)
        self.start_time = 0.
        self.time = 0.
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
