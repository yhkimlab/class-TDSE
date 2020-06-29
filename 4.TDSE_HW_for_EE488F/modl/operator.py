import numpy as np
import numpy.linalg as lin
from modl import operator as oprt
from modl import Input

# Inputs & Unit conversion
k = (Input.ewave*2/27.211)**0.5
pot_shape = Input.pot_shape
pot_height_eV  = Input.pot_height_eV
barrier_thickness = Input.barrier_thickness
algorithm = Input.algorithm
nstep = Input.nstep

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
class wave(grid):
    def __init__(self, L, n, l, dt):
        grid.__init__(self, L, n)
        self.grd = np.array(self.grd, dtype = np.complex64)
        self.t_dt = oprt.time_operator(L, n, l, dt, self.dx) # make a time operator
        self.t_backdt = oprt.time_operatorback(L,n,l,dt,self.dx)
        self.integral = 0.
        self.dx=np.float64(L/n)

        if pot_shape == 4:
            for i in range(n):                                   # wave packet initialize (gaussian)
                if (i > (n*4)//10 and i < (n*6)//10):
                    self.grd[i] = np.exp(-(i*self.dx-0.5*n*self.dx)**2/10)
                else:
                    self.grd[i] = 0. + 0.j
            self.grd /= lin.norm(self.grd)                       # Fix to normalize
            for i in range(n):
                self.grd[i] = self.grd[i]*np.exp(1j*k*(i*self.dx-0.5*n*self.dx))  #Wave packet
        else:
            for i in range(n):                                   # wave packet initialize (gaussian)
                if (i > n*0/10 and i < n*8/10):
                    self.grd[i] = np.exp(-(i*self.dx-0.7*n*self.dx)**2/10)
                else:
                    self.grd[i] = 0. + 0.j
            self.grd /= lin.norm(self.grd)                       # Fix to normalize
            for i in range(n):
                self.grd[i] = self.grd[i]*np.exp(1j*k*(i*self.dx-0.7*n*self.dx))  #Wave packet

# Time propagating
    def next_step(self):
        self.grd = np.dot(self.t_dt.oprt, self.grd)

    def back_step(self):
        self.grd = np.dot(self.t_backdt.oprt, self.grd)

# Construct FDM coefficient
# Define FDM points & coefficients (Here, 7-points FDM)
def fdmcoefficient(l):
    a=np.zeros((2*l+1,2*l+1))
    b=np.zeros(2*l+1)
    c=np.zeros(2*l+1)

    for i in range(0, 2*l+1):
        for j in range(0, 2*l+1):
            a[i,j]= (j-l)**i
    c[2]=2
    a = lin.inv(a)
    b= np.matmul(a, c)

    C=np.zeros((l+1))

    for i in range(0, l+1):
        C[i]=b[i+l]
    return C

# Mold for inheritance.
class Operator:
    def __init__(self, n):
        self.oprt = np.zeros((n, n))
        self.oprt = np.array(self.oprt, dtype = np.float64)

    def oprtprint(self):
        print (self)                       # print Name
        print (self.oprt)                  # print Matrix
        print (self.oprt.dtype)            # print Type

# Construct Laplacian using FDM Coefficient
class Laplacian(Operator):
    def __init__(self, n, l, dx):
        Operator.__init__(self, n)
        self.C = fdmcoefficient(l)    # FDM coefficient
        for i in range(n):
            for j in range(-l, l + 1, 1):
                k = i + j
                if (k >= n):
                    k -= n
                    self.oprt[i][k] = self.C[abs(j)] / (dx**2)
                elif (k<0):
                    k += n
                    self.oprt[i][k] = self.C[abs(j)] / (dx**2)
                else:
                    self.oprt[i][k] = self.C[abs(j)] / (dx**2)


# Construct Potential Operator using Potential grid
class Potential(grid):
    def __init__(self, L, n):
        grid.__init__(self, L, n)
     
        self.grd[0] = 0
        self.grd[1] = 0
        self.grd[n-1]= 0
        #self.grd[0] = 100000000
        #self.grd[1] = 100000000
        #self.grd[n-1]= 100000000
        if pot_shape == 0:                                   #no potential
           self.left = 0.5*n
           self.right = 0.5*n
           for i in range(1, n):
               self.grd[i] = 0                           # Make potential

        if pot_shape == 1:                                   #Step potential
           self.left = 0.5*n
           self.right = 0.5*n
           for i in range(2, n-2):
               self.grd[i] = 0                           # Make potential
           for i in range((n//2),(n-2)):
               self.grd[i] = pot_height_eV              # eV unit
               self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot_shape == 2:                                   #Potential barrier
           self.left = 0.5*n
           self.right = 0.5*n
           for i in range(2, n-2):
               self.grd[i] = 0                           # Make potential
           for i in range((5*n)//10,(5*n)//10+barrier_thickness*10):
               self.grd[i] = pot_height_eV                    # eV unit
               self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot_shape == 3:                                   #Double barrier
            self.left = (45*n)//100-barrier_thickness*10
            self.right = (50*n)//100+barrier_thickness*10
            for i in range(2, n-2):
                self.grd[i] = 0                           # Make potential
            for i in range((45*n)//100-barrier_thickness*10,(45*n)//100):
                self.grd[i] = pot_height_eV                    # eV unit
                self.grd[i] = self.grd[i]/27.211          # eV -> Har
            for i in range((50*n)//100,(50*n)//100+barrier_thickness*10):
                self.grd[i] = pot_height_eV               # eV unit
                self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot_shape == 4:                                   # Square well
            self.left = (n*4)//10
            self.right = (n*6)//10
            for i in range(2, n-2):
                self.grd[i] = 0                           # Make potential
            for i in range(2,(n*4)//10):
                self.grd[i]= pot_height_eV
                self.grd[i] = self.grd[i]/27.211          # eV -> Har
            for i in range((n*6)//10,n-2):
                self.grd[i]= pot_height_eV
                self.grd[i] = self.grd[i]/27.211          # eV -> Har

        if pot_shape == 5:                              #Alpha
            self.left = n//2
            self.right = n//2
            for i in range(0, 500):
                self.grd[i]=0
            for i in range(500, 600):
                self.grd[i]= 66/27.211
            for i in range(600, n):
                self.grd[i]= 50/27.211


        self.oprt = np.zeros((n,n))
        for i in range(0, n):
            self.oprt[i, i]=self.grd[i]
 


# Construct Hamiltonian using Potential and Laplacian
class Hamiltonian(Operator):
    def __init__(self, L, n, l, dx):
        Operator.__init__(self, n)
        self.V = Potential(L, n)
        self.L = Laplacian(n, l, dx)                   # h- = 1, m_e = 1
        self.oprt = -self.L.oprt / 2. + self.V.oprt    # H = - (h-^2/2m) L + V

# Construct Time Operator using Hamiltonian
class time_operator(Operator):
    def __init__(self, L, n, l, dt, dx):
        Operator.__init__(self, n)
        self.oprt        = np.array(self.oprt, dtype = np.complex64)
        self.H           = Hamiltonian(L, n, l, dx)

        if algorithm == 0:                                 #forward method
            self.exp_iHdt_0  = np.eye(n) - 1.j * self.H.oprt * dt
            self.oprt        = np.array(self.exp_iHdt_0, dtype = np.complex64)

        if algorithm == 1:                                 #backward method
            self.exp_iHdt_1  = np.eye(n) + 1.j * self.H.oprt * dt
            self.exp_iHdt_1  = np.array(self.exp_iHdt_1, dtype = np.complex64)
            self.oprt        = lin.inv(self.exp_iHdt_1)
 
        if algorithm == 2:                                 #Crank-Nickolson method
            self.exp_iHdt_2  = np.eye(n) + 1.j * self.H.oprt * (dt) / 2.
            self.exp_iHdt_2  = np.array(self.exp_iHdt_2, dtype = np.complex64)
            self.exp_iHdt_2  = lin.inv(self.exp_iHdt_2)
            self.exp_miHdt_2 = np.zeros_like(self.oprt)
            self.exp_miHdt_2 = np.eye(n) - 1.j * self.H.oprt * (dt) / 2.
            self.oprt        = np.dot(self.exp_iHdt_2, self.exp_miHdt_2)

        print ('dx=', dx)
        print ('dt=', dt)


class time_operatorback(Operator):
    def __init__(self, L, n, l, dt, dx):
        Operator.__init__(self, n)
        self.oprt        = np.array(self.oprt, dtype = np.complex64)
        self.H           = Hamiltonian(L, n, l, dx)
        self.exp_iHdt_0  = np.eye(n) - 1.j * self.H.oprt * dt
        self.oprt        = np.array(self.exp_iHdt_0, dtype = np.complex64)
        self.oprt        = lin.inv(self.oprt)


#               PSI(x,t) = e^-iHt PSI(x,0)
#   e^iHdt/2 PSI(x,t+dt) = e^-iHdt/2 PSI(x,t)
# (1+iHdt/2) PSI(x,t+dt) = (1-iHdt/2) PSI(x,t)
#            PSI(x,t+dt) = (1+iHdt/2)^-1 (a-iHdt/2) PSI(x,t)


