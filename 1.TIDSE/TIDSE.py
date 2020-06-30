import numpy as np
import numpy.linalg as lin
from matplotlib import pyplot as plt
from modl import operator, Input

#Inputs & Unit conversion
pot_height_har = Input.pot_height_eV/27.211
ngx=Input.ngx
l=3
Lx=Input.Lx*1.88973
dx=Lx/(ngx)
nstate = Input.nstate

# Define Hamiltonian using Potential and Laplacian
def Hamiltonian(L,n,l,dx):
    Hamiltonian = np.zeros((n,n))
    V = operator.Potential(L, n)
    L = operator.Laplacian(n, l, dx)                   # h- = 1, m_e = 1
    Hamiltonian = -L.oprt / 2. + V.oprt    # H = - (h-^2/2m) L + V
    return Hamiltonian

# Construct Hamiltonian using Potential and Laplacian
H=np.zeros((ngx, ngx))
dx=Lx/(ngx-1)
dn=np.ones(ngx)/dx**2
dno=np.ones(ngx-1)*(-1/2)/dx**2

V = operator.Potential(Lx,ngx)
H = Hamiltonian(Lx, ngx, l, dx)

#Diagonalize Hamiltonian to get Eigenvalues
(w,v) = np.linalg.eigh(H)
num = np.argsort(w)

#Make box & Visualize
a = np.linspace(-Lx/2, Lx/2, ngx)
a = a/1.89

#Making Total figure
for i in range(0,nstate):
    bb=pot_height_har*27.211/20
    bbb=np.max(v[:,num[i]])-np.min(v[:,num[i]])
    
    plt.plot(a,v[:,num[i]]/bbb*bb+np.ones(ngx)*w[i]*27.211)

plt.plot(a,V.grd*27.211)
plt.ylim(0,Input.pot_height_eV)
plt.ylabel('Energy [eV]')
plt.xlabel('Box [Angstrom]')
plt.savefig('Total.png')
plt.show()

#Making figure for individual eigenstate
for i in range(0,nstate):
    print ('%dth state eigenvalue is %f eV' %(i+1,w[num[i]]*27.211))
    plt.clf()
    plt.plot( a,v[:,num[i]], label='%dth eigenvalue = %f eV' %(i+1,w[num[i]]*27.211))
    plt.legend()
    plt.ylim((-np.max(np.abs(v[:,num[i]])),np.max(np.abs(v[:,num[i]]))))
    plt.savefig('%04d' %(i+1))
