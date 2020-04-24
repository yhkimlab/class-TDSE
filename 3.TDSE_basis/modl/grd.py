import numpy as np
import numpy.linalg as lin
from modl import inp


# Inputs & Unit conversion
k = (inp.E*2/27.211)**0.5
pot = inp.Potential_Shape
pot_height = inp.Potential_Height


# Grid class is real space,, and mold for inheritance.
class grid:
    def __init__(self, L, n):
        self.n   = n
        self.dx  = L/n
        self.grd = np.linspace(-L/2., L/2., num = n, endpoint=False) # Make grid
        self.grd = np.array(self.grd, dtype = np.float64)

# Return string grid data (probability density for waveFunction)
    def plot_grid(self):
        self.strline = ""
        if "float" in str(self.grd.dtype):
            self.tmp_grd = self.grd
        if "complex" in str(self.grd.dtype):
            self.tmp_grd = self.grd.real**2 + self.grd.imag**2
        for i in range(self.n):
            self.strline += "  %.5e " % self.tmp_grd[i]
        self.strline += "  %.6e  \n" % self.tmp_grd[0]
        return self.strline

# wave function : complex grid, construct a time operator
class waveF(grid):
    def __init__(self, L, n, l):
        grid.__init__(self, L, n)
        self.grd = np.array(self.grd, dtype = np.complex64)
        self.integral = 0.
        self.dxx=np.float64(L/n)
        if pot == 3:
            for i in range(n):                                   # wave packet initialize (gaussian)
                if (i > n*0/10 and i < n*9/10):
                    self.grd[i] = np.exp(-(i*self.dxx-0.3*n*self.dxx)**2/9e5)
                else:
                    self.grd[i] = 0. + 0.j
            self.grd /= lin.norm(self.grd)                       # Fix to normalize
            for i in range(n):
                self.grd[i] = self.grd[i]*np.exp(1j*k*(i*self.dxx-0.3*n*self.dxx))  #Wave packet

        if pot == 4:
            for i in range(n):                                   # wave packet initialize (gaussian)
                if (i > n*4/10 and i < n*6/10):
                    self.grd[i] = np.exp(-(i*self.dxx-0.5*n*self.dxx)**2/10)
                else:
                    self.grd[i] = 0. + 0.j
            self.grd /= lin.norm(self.grd)                       # Fix to normalize
            for i in range(n):
                self.grd[i] = self.grd[i]*np.exp(1j*k*(i*self.dxx-0.5*n*self.dxx))  #Wave packet

        else:
            for i in range(n):                                   # wave packet initialize (gaussian)
                if (i > n*0/10 and i < n*4/10):
                    self.grd[i] = np.exp(-(i*self.dxx-0.3*n*self.dxx)**2/10)
                else:
                    self.grd[i] = 0. + 0.j
            self.grd /= lin.norm(self.grd)                       # Fix to normalize
            for i in range(n):
                self.grd[i] = self.grd[i]*np.exp(1j*k*(i*self.dxx-0.3*n*self.dxx))  #Wave packet

# Set the potential
class Potential(grid):
    def plot_grid(self):
        self.strline = ""
        if "float" in str(self.grd.dtype):
            self.tmp_grd = self.grd
        if "complex" in str(self.grd.dtype):
            self.tmp_grd = self.grd.real**2 + self.grd.imag**2
        for i in range(self.n):
            self.strline += "  %.5e " % self.tmp_grd[i]
        self.strline += "  %.6e  \n" % self.tmp_grd[0]
        return self.strline

    def __init__(self, L, n):
        grid.__init__(self, L, n)
        self.grd[0] = 100000000
        self.grd[1] = 100000000
        self.grd[n-1]= 100000000
        if pot == 0:                                   #no potential
           self.left = 0.5*n
           self.right = 0.5*n
           for i in range(1, n):
               self.grd[i] = 0                           # Make potential

        if pot == 1:                                   #Step potential
           self.left = 0.5*n
           self.right = 0.5*n
           for i in range(2, n-2):
               self.grd[i] = 0                           # Make potential
           for i in range((n//2),(n-2)):
               self.grd[i] = pot_height                   # eV unit
               self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot == 2:                                   #Potential barrier
           self.left = 0.5*n
           self.right = 0.5*n
           for i in range(2, n-2):
               self.grd[i] = 0                           # Make potential
           for i in range(450,551):
               self.grd[i] = pot_height                   # eV unit
               self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot == 3:                                   #Double barrier
            self.left = 0.6*n
            self.right = n*0.6
            for i in range(2, n-2):
                self.grd[i] = 0                           # Make potential
            for i in range(1200, 1204):
                self.grd[i] = pot_height                   # eV unit
                self.grd[i] = self.grd[i]/27.211          # eV -> Har
            for i in range(1209,1213):
                self.grd[i] = pot_height                   # eV unit
                self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot == 4:                                   # Square well
            self.left = 0.5*n
            self.right = 0.5*n
            for i in range(2, n-2):
                self.grd[i] = 0                           # Make potential
            for i in range(2,(n*4)//10):
                self.grd[i]= pot_height
                self.grd[i] = self.grd[i]/27.211          # eV -> Har
            for i in range((n*6)//10,n-2):
                self.grd[i]= pot_height
                self.grd[i] = self.grd[i]/27.211          # eV -> Har
