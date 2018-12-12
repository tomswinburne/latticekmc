import numpy as np
import matplotlib.pyplot as plt
import time,sys
sys.path.insert(0,'./lib/')
from engine import KMCEngine

Etheta = np.pi/3.
qE = np.r_[[np.cos(Etheta),np.sin(Etheta)]] * 6. # in kT

kmc = KMCEngine(lattice='honeycomb',jump=2.0,bond=2.,force=qE,penalty=0.)

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
if len(sys.argv)>2:
    temp = np.loadtxt(sys.argv[2]).astype(bool)
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
if len(sys.argv)>1:
    plt.savefig(sys.argv[1])
else:
    plt.show()
