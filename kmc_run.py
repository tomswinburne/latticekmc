import numpy as np
import matplotlib.pyplot as plt
import time,sys
sys.path.insert(0,'./c++/')
from engine import KMCEngine

simtype = 'traj'


jump = 1.0
bond = 3.0
penalty = 0.
theta=0.#+np.pi/8.
#force = np.r_[[np.cos(theta),np.sin(theta)]]*2.

force = np.r_[[2.,0.]] + np.ones(2)


"""
basis = np.zeros((4,2))
cell = np.r_[[2.*np.cos(np.pi/3.)+2.,2.*np.sin(np.pi/3.)]]
basis[1][0] = 1.
basis[1][1] = 0.
basis[2][0] = np.cos(np.pi/3.)+1.
basis[2][1] = np.sin(np.pi/3.)
basis[3][0] = np.cos(np.pi/3.)+2.
basis[3][1] = np.sin(np.pi/3.)
rad = 1.1#*np.sqrt( 2. + 2.*np.cos(np.pi/3.) ) * 1.1
"""


"""
cell = np.r_[[1.,2.*np.sin(np.pi/3.)]]
basis=np.zeros((2,2))
basis[1][0] = np.cos(np.pi/3.)
basis[1][1] = np.sin(np.pi/3.)
rad = 1.1
"""

cell = np.r_[[1.,1.]]#2.*np.sin(np.pi/3.)]]
basis=np.zeros((2,2))
basis[1][0] = .5#np.cos(np.pi/4.)
basis[1][1] = .5#np.sin(np.pi/4.)
rad = 0.9



ncells = [200,100]
box = [int(ncells[0]/cell[0])*cell[0],int(ncells[1]/cell[1])*cell[1]]
cc = np.r_[[box[0]-4.*rad,box[1]/2.-4.*rad]]
slopes = [np.tan(-np.pi/8.),np.tan(np.pi/4.),-np.tan(np.pi/4.)]




sim = KMCEngine(cell=cell, basis=basis, radius=rad,\
                jump=jump, bond=bond, force=force, penalty=penalty)
sim.build(ncells[0],ncells[1])

pos = sim.get_positions()
center =  pos[np.linalg.norm(pos-cc,axis=1).argmin()]
sel = np.linalg.norm(pos-center,axis=1)>rad*4
if len(sys.argv) > 1:
    occupation = np.loadtxt(sys.argv[1]).astype(bool).flatten()
else:
    occupation = np.zeros(sim.nsites,bool)
    occupation[sel] = True

ftraj = None
if len(sys.argv) > 2:
    ftraj = np.loadtxt(sys.argv[2])
    sim.start_time = ftraj[-1][0]

sim.set_occupations(occupation)

def line(slope,point):
    l = np.zeros((100,2))
    l[:,0] = np.linspace(0.,box[0],100)
    l[:,1] = slope*(l[:,0]-point[0])+point[1]
    return l

c = np.zeros((2,2))
ic = np.zeros((2,2))
for i in range(2):
    c[i][i] = box[i]
    ic[i][i] = 1./box[i]
def pbc(_pos):
    # find dominant quadrant and pick point in quad
    ca = 0
    _pbc = (_pos-_pos[ca]-np.round((_pos-_pos[ca]).dot(ic)).dot(c))
    _pbcd = np.linalg.norm(_pbc,axis=1)
    while np.percentile(_pbcd,25) > 6.*rad:
        ca += 1
        _pbc = (_pos-_pos[ca]-np.round((_pos-_pos[ca]).dot(ic)).dot(c))
        _pbcd = np.linalg.norm(_pbc,axis=1)
        if ca==_pos.shape[0]-1:
            print "END"
            break
    return _pbc[_pbcd<6.*rad].mean(axis=0)+_pos[ca]



if simtype == "zoetrope":
    substeps = 40
    nframes = 2
    yf = 1
    xf = nframes/yf
    nframes = xf*yf
    fig = plt.figure(figsize=(2*xf,2*yf*box[1]/box[0]))
    traj = []
    for sp in range(nframes):
        ax = plt.subplot(yf,xf,sp+1)
        if sp>0:
            for ss in range(substeps):
                sim.run(50000,True,sp*substeps+ss)
                occ = sim.get_occupations().flatten()
                traj.append(pbc(pos[~occ]))
        else:
            ll = 3./np.linalg.norm(force)
            ax.arrow(3,box[1]/2.,force[0]*ll,force[1]*ll,width=1.)
            occ = sim.get_occupations().flatten()
            traj.append(pbc(pos[~occ]))



        ax.set_xlim(0,box[0])
        ax.set_ylim(0,box[1])
        #occ = sim.get_occupations().flatten()
        nnv = sim.get_neigh_count().flatten()
        ax.scatter(pos[:,0], pos[:,1],s=3,c=nnv,vmin=0)

        #line += center
        l = line(np.tan(-np.pi/3.),center)
        ax.plot(l[:,0],l[:,1],'k--',lw=1.)
        l = line(np.tan(-np.pi/6.),center)
        ax.plot(l[:,0],l[:,1],'k--',lw=1.)
        l = line(np.tan(np.pi/3.),center)
        ax.plot(l[:,0],l[:,1],'k--',lw=1.)
        if sp>0:
            _traj = np.r_[traj]
            ax.scatter(_traj[:,0],_traj[:,1],c='k',s=4)

    np.savetxt('occupation',sim.get_occupations().flatten())
    plt.show()


if simtype == "traj":
    inspections = 4000
    steps = 1000
    _traj = np.zeros((inspections,3))
    occ = sim.get_occupations().flatten()
    _traj[0][0] = sim.get_time()
    _traj[0][1:3] = pbc(pos[~occ])
    now = time.time()
    for i in range(1,inspections):
        sim.run(steps=steps,restart=False,seed=i)
        occ = sim.get_occupations().flatten()
        _traj[i][0] = sim.get_time()
        _traj[i][1:3] = pbc(pos[~occ])
        if i%100==0:
            print "Ran for %d steps in %2.2fs, %d/%d" % \
                (100*steps,time.time()-now,i,inspections)
            now = time.time()

    fig = plt.figure(figsize=(6,6*box[1]/box[0]))
    ax = plt.subplot(1,1,1)
    ll = 3./np.linalg.norm(force)
    ax.arrow(3,box[1]/2.,force[0]*ll,force[1]*ll,width=1.)
    ax.set_xlim(0,box[0])
    ax.set_ylim(0,box[1])


    nnv = sim.get_neigh_count().flatten()
    ax.scatter(pos[:,0], pos[:,1],s=5,c=nnv,vmin=0)
    for s in slopes:
        l = line(s,center)
        ax.plot(l[:,0],l[:,1],'k--',lw=1.)
    if not ftraj is None:
        traj = np.vstack((ftraj,_traj))
    else:
        traj = _traj.copy()
    ax.plot(traj[:,1],traj[:,2],'k-')#c='k',s=4)
    np.savetxt('occupation',sim.get_occupations().flatten())
    np.savetxt('traj',traj)
    plt.show()



exit()




#print pos.shape
#sim.run(10000)
#exit()

#print pos.shape,occupation.shape
#thet = np.linspace(0.,2.*np.pi,100)
occupation = sim.get_occupations()
fig = plt.figure(figsize=(12,6))
ax = plt.subplot(1,2,1)
pos = sim.get_positions()[occupation.flatten()]
ax.scatter(pos[:,0],pos[:,1],s=1.)
sim.run(100000)

ax2 = plt.subplot(1,2,2)
occupation = sim.get_occupations()
pos = sim.get_positions()[occupation.flatten()]
ax2.scatter(pos[:,0],pos[:,1],s=1.)
plt.tight_layout()
plt.show()



#sim.test(5.0)


"""

oo = np.linspace(0.,1.,10)
oob = oo.astype(np.bool)
print oo[0]






Etheta = 0.#np.pi/3.
qE = np.r_[[np.cos(Etheta),np.sin(Etheta)]] * 1. # in kT

kmc = KMCEngine(lattice='honeycomb',jump=1.0,bond=3.,force=qE,penalty=10.)

kmc.build(xsize=50,ysize=50)
ratio = 1.
nframes = 12
yf = 3

# cut out a circle - put this in a function?
center = kmc.positions.mean(axis=0)
center = kmc.positions[np.linalg.norm(kmc.positions - center,axis=1).argmin()]
sel = np.linalg.norm(kmc.positions - center,axis=1) < 5.
kmc.occupation[sel] = False


# load in from file (optional)
if len(sys.argv)>1:
    temp = np.loadtxt(sys.argv[1]).astype(bool)
    kmc.occupation = temp.copy()
else:
    np.savetxt('start',kmc.occupation)


xf = nframes/yf
nframes = xf*yf
fig = plt.figure(figsize=(3*xf,3*yf*ratio))
for sp in range(nframes):
    ax = plt.subplot(yf,xf,sp+1)
    if sp>0:
        kmc.run(steps=5000, restart=False)
    else:
        ll = 3. / np.linalg.norm(qE)
        ax.arrow(3,kmc.height-3,qE[0]*ll,qE[1]*ll,width=1.)

    print "%10f %d/%d" % (kmc.time,sp+1,nframes)
    #ax.set_title("time = %4.4g" % sim_time)
    ax.set_xlim(0,kmc.Ncells[0]*kmc.cell[0])
    ax.set_ylim(0,kmc.Ncells[1]*kmc.cell[1])
    num_neigh = kmc.occupation[kmc.neighbors].sum(axis=1)[kmc.occupation]
    rpos = kmc.positions[kmc.occupation]
    ax.scatter(rpos[:,0], rpos[:,1],\
                    s=5,c=num_neigh,cmap="viridis",vmin=0,vmax=kmc.moves)

np.savetxt('occupation',kmc.occupation)

plt.tight_layout()
if len(sys.argv)>2:
    plt.savefig(sys.argv[2])
else:
    plt.show()

"""
