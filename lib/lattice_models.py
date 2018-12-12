import numpy as np

"""
Current restriction is to deliver orthorhombic unit cell to code
This lattice class could be given biclinic cell then generate orthorhombic...
"""
class lattice_model:
    def __init__(self,motif="hcp"):
        """ Define primitive cell (potentially biclinic in future) """
        self.primitive_cell = np.identity(2)

        if motif == "hcp":
            self.primitive_cell[1][1] = 2.0*np.sin(np.pi/3.)
            self.basis=np.zeros((2,2))
            self.basis[1][0] = np.cos(np.pi/3.)
            self.basis[1][1] = np.sin(np.pi/3.)
            self.rad = 1.1

        elif motif == "fcc":
            self.basis=np.zeros((2,2))
            self.basis[1][0] = .5
            self.basis[1][1] = .5
            self.rad = 1.1/np.sqrt(2.)

        elif motif == "graphene":
            self.primitive_cell[0][0] = 2.0*(1+np.cos(np.pi/3.))
            self.primitive_cell[1][1] = 2.0*(0+np.sin(np.pi/3.))
            self.basis = np.zeros((4,2))
            self.basis[1][0] = 1.
            self.basis[1][1] = 0.
            self.basis[2][0] = np.cos(np.pi/3.)+1.
            self.basis[2][1] = np.sin(np.pi/3.)
            self.basis[3][0] = np.cos(np.pi/3.)+2.
            self.basis[3][1] = np.sin(np.pi/3.)
            self.rad = 1.1
        else:
            print "Unknown lattice type"
            exit(-1)


    def unit_cell(self):
        # Currently trivial as primitive_cell is orthorhombic
        unit_cell = self.primitive_cell.diagonal()
        return unit_cell

    def unit_basis(self):
        # Currently trivial as primitive_cell is orthorhombic
        unit_basis = self.basis
        return unit_basis
