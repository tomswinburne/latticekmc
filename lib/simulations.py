import sys,time
import numpy as np
import matplotlib.pyplot as plt




class simulator:
    def __init__(self, sim, cycles=1000, steps_per_cycle=1000,\
                snapshots=10, file=None, lines=None,\
                occupation_file=None,traj_file=None, cluster_file=None):
        self.sim = sim
        self.cycles = cycles
        self.steps_per_cycle = steps_per_cycle

        self.file = file
        self.traj_file = traj_file
        self.cluster_file = cluster_file
        self.occupation_file = occupation_file
        self.lines = lines
        self.snapshots = snapshots

        """ lattice positions used throughout """
        self.pos = sim.get_positions()


        """ overwrite if matching file exists """
        if not self.occupation_file is None:
            try:
                test_occ = np.loadtxt(occupation_file).astype(bool).flatten()
                if test_occ.shape[0] == self.sim.nsites:
                    self.sim.set_occupations(test_occ)
                    print "Loaded occupation file"
            except IOError:
                print "Can't find occupation file"

        """ cluster position function """
        self.c = np.zeros((2,2))
        self.ic = np.zeros((2,2))
        for i in range(2):
            self.c[i][i] = sim.super_cell[i]
            self.ic[i][i] = 1./sim.super_cell[i]

        """ load in trajectory file, if it exists """
        self.ftraj = None
        if not self.traj_file is None:
            try:
                self.ftraj = np.loadtxt(traj_file)
                sim.start_time = self.ftraj[-1][0]
                print "Read in traj file, start time:",self.ftraj[-1][0]
            except IOError:
                self.traj_file = None

        """ load in cluster file, if it exists """
        self.ctraj = None
        if not self.cluster_file is None:
            try:
                self.ctraj = np.loadtxt(self.cluster_file)
                print "Read in cluster file"
            except IOError:
                self.cluster_file = None

    def cluster_position(self,_pos):
        ca = 0
        _pbc = _pos-_pos[ca]
        _pbc -= np.round(_pbc.dot(self.ic)).dot(self.c)
        _pbcd = np.linalg.norm(_pbc,axis=1)

        while np.percentile(_pbcd,25) > 6.*self.sim.radius:
            ca += 1
            _pbc = _pos-_pos[ca]
            _pbc -= np.round(_pbc.dot(self.ic)).dot(self.c)
            _pbcd = np.linalg.norm(_pbc,axis=1)
            if ca==_pos.shape[0]-1:
                break
        return _pbc[_pbcd<6.*self.sim.radius].mean(axis=0)+_pos[ca]

    def run(self,type='traj'):
        """ trajectory for this run """
        traj = np.zeros((self.cycles,3))
        occ = self.sim.get_occupations().flatten()
        traj[0][0] = self.sim.get_time()
        traj[0][1:3] = self.cluster_position(self.pos[~occ])

        cluster_nnv = []
        now = time.time()
        for cyc in range(1,self.cycles):

            """ run for cycles """
            self.sim.run(steps=self.steps_per_cycle,restart=False,seed=cyc)

            """ add to trajectory """
            occ = self.sim.get_occupations().flatten()
            traj[cyc][0] = self.sim.get_time()
            traj[cyc][1:3] = self.cluster_position(self.pos[~occ])

            """ print to screen """
            if cyc%100==99:
                print "Ran for %d steps in %2.2fs, %d/%d" % \
                    (100*self.steps_per_cycle,time.time()-now,cyc+1,self.cycles)
                now = time.time()

            if cyc % (self.cycles/self.snapshots)==0:
                print "SNAPSHOT"
                nnv = self.sim.get_neigh_count().flatten()[~occ].reshape((-1,1))
                cluster_nnv.append(np.hstack((self.pos[~occ],nnv)))


        """ plot snapshots and traj """
        ratio = self.sim.super_cell[1] / self.sim.super_cell[0]
        fig = plt.figure(figsize=(6,12*ratio))

        ax = plt.subplot(2,1,1)
        ax.set_xlim(0,self.sim.super_cell[0])
        ax.set_ylim(0,self.sim.super_cell[1])

        """ force arrow """
        if np.linalg.norm(self.sim.force)>0.001:
            arrow = 3./np.linalg.norm(self.sim.force) * self.sim.force
            ax.arrow(3,self.sim.super_cell[1]/2.,arrow[0],arrow[1],width=1.)

        """ Guide lines """
        if not self.lines is None:
            for sl_p in self.lines:
                xx = np.linspace(0.,self.sim.super_cell[0],100)
                yy = (xx-sl_p[1][0])*sl_p[0] + sl_p[1][1]
                ax.plot(xx,yy,'k--',lw=1.)

        """ Plot trajectory """
        if self.ftraj is None:
            _traj = traj.copy()
        else:
            _traj = np.vstack((self.ftraj,traj))

        ax.scatter(_traj[:,1],_traj[:,2],c='k',s=4)



        if self.cluster_file is None:
            self.ctraj = cluster_nnv[0]
        else:
            self.ctraj = np.vstack((self.ctraj,cluster_nnv[0]))
        for cv in cluster_nnv[1:]:
            self.ctraj = np.vstack((self.ctraj,cv))

        ax.scatter(self.ctraj[:,0],self.ctraj[:,1],s=5,c=self.ctraj[:,2],vmin=0)

        xy = _traj[:,[1,2]]
        xy -= xy[0]
        dxy = xy[1:]-xy[:-1]
        xy[1:] = (dxy-np.round(dxy.dot(self.ic)).dot(self.c)).cumsum(axis=0)
        ax2 = plt.subplot(2,1,2)
        ax2.set_xlim(xy[:,0].min(),xy[:,0].max())
        ax2.set_ylim(xy[:,1].min(),xy[:,1].max())
        ax2.plot(xy[:,0],xy[:,1],'C1-',lw=2)

        """ Guide lines """
        if not self.lines is None:
            for sl_p in self.lines:
                xx = np.linspace(xy[:,0].min(),xy[:,0].max(),100)
                ax2.plot(xx,xx*sl_p[0],'k--',lw=1.)

        """ Save configuration, trajectory and cluster """
        if not self.occupation_file is None:
            np.savetxt(self.occupation_file,\
                self.sim.get_occupations().flatten())
        else:
            np.savetxt('occupation',self.sim.get_occupations().flatten())

        if not self.traj_file is None:
            np.savetxt(self.traj_file,_traj)
        else:
            np.savetxt('traj',_traj)


        if not self.cluster_file is None:
            np.savetxt(self.cluster_file,self.ctraj)
        else:
            np.savetxt('cluster',self.ctraj)

        """ Show or save plot """
        if self.file is None:
            plt.show()
        else:
            plt.savefig(self.file)
