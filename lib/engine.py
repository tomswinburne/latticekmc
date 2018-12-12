import numpy as np
class KMCEngine:
    def __init__(self,lattice,jump=1.0,bond=4.0,force=np.zeros(2),penalty=1.):
        self.jump = jump
        self.bond = bond
        self.penalty = penalty
        self.force = force
        self.moves=0
        self.nmoves=False
        self.time=0.
        if lattice == 'honeycomb':
            self.moves = 3
            self.nmoves = 6
            self.basis = np.zeros((4,2))
            self.cell = np.r_[[2.*np.cos(np.pi/3.)+2.,2.*np.sin(np.pi/3.)]]

            self.basis[1][0] = 1.
            self.basis[1][1] = 0.
            self.basis[2][0] = np.cos(np.pi/3.)+1.
            self.basis[2][1] = np.sin(np.pi/3.)
            self.basis[3][0] = np.cos(np.pi/3.)+2.
            self.basis[3][1] = np.sin(np.pi/3.)
        elif lattice == 'hcp':
            self.moves = 6
            self.nmoves = 0
            self.basis = np.zeros((2,2))
            self.cell = np.r_[[1.,2.*np.sin(np.pi/3.)]]

            self.basis[1][0] = np.cos(np.pi/3.)
            self.basis[1][1] = np.sin(np.pi/3.)
        else:
            print "Unknown Lattice"

        self.Nbasis = len(self.basis)
        self.nlist = np.zeros((self.Nbasis,self.moves,3),dtype=int)
        if self.nmoves:
            self.nnlist = np.zeros((self.Nbasis,self.nmoves,3),dtype=int)
        self.jvec = np.zeros((self.Nbasis,self.moves,2))

        if lattice=='honeycomb':
            self.nlist[0] = np.r_[[[0,0,1],[-1,0,3],[-1,-1,3]]]
            self.nlist[1] = np.r_[[[0,0,0],[0,0,2],[0,-1,2]]]
            self.nlist[2] = np.r_[[[0,0,1],[0,1,1],[0,0,3]]]
            self.nlist[3] = np.r_[[[0,0,2],[1,0,0],[1,1,0]]]

            self.nnlist[0] = np.r_[[[0,1,0],[0,-1,0],[0,0,2],\
                                    [-1,0,2],[0,-1,2],[-1,-1,2]]]
            self.nnlist[1] = np.r_[[[0,1,1],[0,-1,1],[0,0,3],\
                                    [-1,0,3],[0,-1,3],[-1,-1,3]]]
            self.nnlist[2] = np.r_[[[0,1,2],[0,1,2],[0,0,0],\
                                    [1,0,0],[0,1,0],[1,1,0]]]
            self.nnlist[3] = np.r_[[[0,0,1],[1,0,1],[0,1,1],\
                                    [1,1,1],[0,1,3],[0,-1,3]]]
        elif lattice=='hcp':
            self.nlist[0] = np.r_[[[0,0,1],[1,0,0],[-1,0,0],\
                                [-1,0,1],[0,-1,1],[-1,-1,1]]]
            self.nlist[1] = np.r_[[[0,0,0],[1,0,1],[-1,0,1],\
                                [1,0,0],[0,1,0],[1,1,0]]]

        for s in range(self.Nbasis):
            for i,nl in enumerate(self.nlist[s]):
                for x in range(2):
                    self.jvec[s][i][x] = nl[x]*self.cell[x] + \
                        self.basis[nl[2]][x] - self.basis[s][x]
            # if self.nmoves....


    def build(self,xsize=10,ysize=None):
        if ysize is None:
            ysize = xsize
        self.Ncells = [int(xsize/self.cell[0]),int(ysize/self.cell[1])]
        self.Nsites = self.Ncells[0] * self.Ncells[1] * self.Nbasis
        self.neighbors = np.zeros((self.Nsites,self.moves),dtype=int)
        if self.nmoves:
            self.next_neighbors = np.zeros((self.Nsites,self.nmoves),dtype=int)
        self.nm = np.ones(self.moves)
        self.indicies = np.outer(np.arange(self.Nsites),self.nm).astype(int)

        self.height = self.Ncells[1] * self.cell[1]
        self.width = self.Ncells[0] * self.cell[0]

        # Single O(N) loop to start
        for i in range(self.Nsites):
            bi = i % self.Nbasis
            ci = i / self.Nbasis
            for ni,n in enumerate(self.nlist[bi]):
                cx = (ci%self.Ncells[0]+n[0]+self.Ncells[0])%self.Ncells[0]
                cy = (ci/self.Ncells[0]+n[1]+self.Ncells[1])%self.Ncells[1]
                self.neighbors[i][ni] = (cx+self.Ncells[0]*cy)*self.Nbasis+n[2]
        if self.nmoves:
            for i in range(self.Nsites):
                bi = i % self.Nbasis
                ci = i / self.Nbasis
                for ni,n in enumerate(self.nnlist[bi]):
                    cx = (ci%self.Ncells[0]+n[0]+self.Ncells[0])%self.Ncells[0]
                    cy = (ci/self.Ncells[0]+n[1]+self.Ncells[1])%self.Ncells[1]
                    self.next_neighbors[i][ni] = \
                                (cx + self.Ncells[0]*cy) * self.Nbasis + n[2]

        # make positions, occupation vector etc
        self.type = np.arange(self.Nsites, dtype=int) % self.Nbasis
        self.occupation = np.ones(self.Nsites, dtype=bool)
        pos = np.zeros((self.Nsites/self.Nbasis,self.Nbasis,2))
        cr = np.arange(self.Nsites/self.Nbasis).astype(int)
        for i,b in enumerate(self.basis):
            pos[:,i,0] = self.cell[0] * (cr % self.Ncells[0]) + b[0]
            pos[:,i,1] = self.cell[1] * (cr / self.Ncells[0]) + b[1]

        self.positions = pos.reshape((-1,2))


    def cycle(self):
        # number of neighbors
        nn = self.occupation[self.neighbors].sum(axis=1)
        if self.nmoves:
            nnn = self.occupation[self.neighbors].sum(axis=1)

        # mobile atoms
        mobile = ( nn < self.moves ) * self.occupation
        mobile_index = self.indicies[mobile]
        # neighbors of mobile atoms
        neighbors = self.neighbors[mobile]
        num_neigh = nn[mobile]

        # list of neighbors which are unoccupied
        valid_dest = ~self.occupation[neighbors.flatten()]
        valid_dest_index = neighbors.flatten()[valid_dest]
        valid_origin_index = mobile_index.flatten()[valid_dest]

        # neighboring sites of each valid destination
        nsites_dest = self.neighbors[neighbors.flatten()][valid_dest]
        nsites_orig = self.neighbors[mobile_index.flatten()][valid_dest]
        #print nsites_dest.shape,nsites_dest.shape
        #print np.outer(nsites_dest,nsites_orig).shape
        #print nsites_dest.shape,nsites_orig.shape
        #,nn[nsites].shape

        # jump vectors
        jvec_dest = self.jvec[self.type[mobile]].reshape((-1,2))[valid_dest]

        # num_neigh of destination: -1 accounts for jumping atom
        nn_dest = nn[neighbors].flatten()[valid_dest]-1

        # num_neigh currently
        nn_cur =  np.outer(nn[mobile],self.nm).flatten()[valid_dest]

        # next neighbors
        if self.nmoves:
            nnn_dest = nnn[neighbors].flatten()[valid_dest]-1
            nnn_cur = np.outer(nn[mobile],self.nm).flatten()[valid_dest]

        # i->j gives n_i atoms one less neigh and n_j atoms one more: -
        #   1*n_i+1*n_j
        # i-> gives jumping atom n_j-n_i more neighbors
        penalty = np.zeros(nn_dest.shape)
        #penalty += 4.*(nn_dest<1)*(nn_cur>0)
        #penalty += self.penalty*(nn_dest>0)*(nn_cur==0)
        penalty += 2.*self.penalty*(nn_cur>0)*(nn_dest>0)

        dE_ij = -2.*self.bond*(nn_dest-nn_cur) - jvec_dest.dot(self.force)
        if self.nmoves:
            dE_ij += -self.bond*(nnn_dest-nnn_cur)
        dE_ij += penalty

        barrier = dE_ij*(dE_ij>0) + self.jump #+ np_penalty

        rate_vec = np.exp(-barrier)
        tr = rate_vec.sum()
        rate_vec /= tr
        u = np.random.uniform(1.e-20,1.,size=2)
        self.time += -np.log(u[1]) / tr
        move = ( rate_vec.cumsum() < u[0] ).argmin()

        old_index = valid_origin_index[move]
        new_index = valid_dest_index[move]
        self.occupation[old_index] = False
        self.occupation[new_index] = True

    def run(self,steps=100,restart=False):
        if restart:
            self.time=0.
        for i in range(steps):
            self.cycle()
