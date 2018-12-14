import sys,time
import numpy as np
import matplotlib.pyplot as plt

def run_and_plot(sim,cycles=1000,steps_per_cycle=1000,snapshots=10,lines=None,\
        image_file=None,occupation_file=None,traj_file=None,cluster_file=None):

        """ lattice positions used throughout """
        pos = sim.get_positions()


        """ overwrite if matching file exists """
        if not occupation_file is None:
            try:
                test_occ = np.loadtxt(occupation_file).astype(bool).flatten()
                if test_occ.shape[0] == sim.nsites:
                    sim.set_occupations(test_occ)
                    print "Loaded occupation file"
            except IOError:
                print "No occupation file found"

        """ cluster position function """
        c = np.zeros((2,2))
        ic = np.zeros((2,2))
        for i in range(2):
            c[i][i] = sim.super_cell[i]
            ic[i][i] = 1./sim.super_cell[i]
        def cluster_position(_pos):
            ca = 0
            _pbc = _pos-_pos[ca]
            _pbc -= np.round(_pbc.dot(ic)).dot(c)
            _pbcd = np.linalg.norm(_pbc,axis=1)

            while np.percentile(_pbcd,25) > 6.*sim.radius:
                ca += 1
                _pbc = _pos-_pos[ca]
                _pbc -= np.round(_pbc.dot(ic)).dot(c)
                _pbcd = np.linalg.norm(_pbc,axis=1)
                if ca==_pos.shape[0]-1:
                    break
            return _pbc[_pbcd<6.*sim.radius].mean(axis=0)+_pos[ca]

        """ load in trajectory file, if it exists """
        ftraj = None
        if not traj_file is None:
            try:
                ftraj = np.loadtxt(traj_file)
                sim.start_time = ftraj[-1][0]
                print "Read in traj file, start time:",ftraj[-1][0]
            except IOError:
                print "No traj file found"

        """ load in cluster file, if it exists """
        ctraj = None
        if not cluster_file is None:
            try:
                ctraj = np.loadtxt(cluster_file)
                print "Read in cluster file"
            except IOError:
                print "No cluster file found"

        """ trajectory for this run """
        traj = np.zeros((cycles,3))
        occ = sim.get_occupations().flatten()
        traj[0][0] = sim.get_time()
        traj[0][1:3] = cluster_position(pos[~occ])

        cluster_nnv = []
        now = time.time()
        for cyc in range(1,cycles):

            """ run for cycles """
            sim.run(steps=steps_per_cycle,restart=False,seed=cyc)

            """ add to trajectory """
            occ = sim.get_occupations().flatten()
            traj[cyc][0] = sim.get_time()
            traj[cyc][1:3] = cluster_position(pos[~occ])

            """ print to screen """
            if cyc%100==99:
                print "Ran for %d steps in %2.2fs, %d/%d" % \
                    (100*steps_per_cycle,time.time()-now,cyc+1,cycles)
                now = time.time()

            if cyc % (cycles/snapshots)==0:
                print "Made Snapshot"
                nnv = sim.get_neigh_count().flatten()[~occ].reshape((-1,1))
                cluster_nnv.append(np.hstack((pos[~occ],nnv)))


        """ plot snapshots and traj """
        ratio = sim.super_cell[1] / sim.super_cell[0]
        fig = plt.figure(figsize=(6,12*ratio))

        ax = plt.subplot(2,1,1)
        ax.set_xlim(0,sim.super_cell[0])
        ax.set_ylim(0,sim.super_cell[1])

        """ force arrow """
        if np.linalg.norm(sim.force)>0.001:
            arrow = 3./np.linalg.norm(sim.force) * sim.force
            ax.arrow(3,sim.super_cell[1]/2.,arrow[0],arrow[1],width=1.)

        """ Guide lines """
        if not lines is None:
            for sl_p in lines:
                xx = np.linspace(0.,sim.super_cell[0],100)
                yy = (xx-sl_p[1][0])*sl_p[0] + sl_p[1][1]
                ax.plot(xx,yy,'k--',lw=1.)

        """ Plot trajectory """
        if ftraj is None:
            _traj = traj.copy()
        else:
            _traj = np.vstack((ftraj,traj))

        ax.scatter(_traj[:,1],_traj[:,2],c='k',s=4)

        if ctraj is None:
            ctraj = cluster_nnv[0]
        else:
            ctraj = np.vstack((ctraj,cluster_nnv[0]))
        for cv in cluster_nnv[1:]:
            ctraj = np.vstack((ctraj,cv))

        ax.scatter(ctraj[:,0],ctraj[:,1],s=5,c=ctraj[:,2],vmin=0)

        xy = _traj[:,[1,2]]
        xy -= xy[0]
        dxy = xy[1:]-xy[:-1]
        xy[1:] = (dxy-np.round(dxy.dot(ic)).dot(c)).cumsum(axis=0)
        ax2 = plt.subplot(2,1,2)
        ax2.set_xlim(xy[:,0].min(),xy[:,0].max())
        ax2.set_ylim(xy[:,1].min(),xy[:,1].max())
        ax2.plot(xy[:,0],xy[:,1],'C1-',lw=2)

        """ Guide lines """
        if not lines is None:
            for sl_p in lines:
                xx = np.linspace(xy[:,0].min(),xy[:,0].max(),100)
                ax2.plot(xx,xx*sl_p[0],'k--',lw=1.)

        """ Save configuration, trajectory and cluster """
        if not occupation_file is None:
            np.savetxt(occupation_file,\
                sim.get_occupations().flatten())
        else:
            np.savetxt('occupation',sim.get_occupations().flatten())

        if not traj_file is None:
            np.savetxt(traj_file,_traj)
        else:
            np.savetxt('traj',_traj)


        if not cluster_file is None:
            np.savetxt(cluster_file,ctraj)
        else:
            np.savetxt('cluster',ctraj)

        """ Show or save plot """
        if image_file is None:
            plt.show()
        else:
            plt.savefig(image_file)
