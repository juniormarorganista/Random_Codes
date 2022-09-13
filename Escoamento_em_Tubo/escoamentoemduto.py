import numpy as np
import math
import time
import os
import shutil
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import diags
from mshr import *
from numpy import linalg
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import inv
from uvw import RectilinearGrid, DataArray

#IO: Creating dir
cwd =os.getcwd()
dir_=cwd+'/Theta0.5-mu1e-3-DT1e-8-33-33-T1e-3'

if not os.path.exists(dir_):
    os.makedirs(dir_)

print(dir_)

### Parameters time
#(0, 0.5 e 1) 
theta = 0.5
tn    = 0.0
tf    = 1.e-4
dt    = 1e-8
T     = 1e-3
nt    = int((tf-tn)/dt)

### Parameters
ro    = 1000
mu    = 1e-3
Lw    = 20.e-6
Lh    = 10.e-6
Lz    = 3.e-3
N1    = 33
N2    = 33
N3    = 33
h1    = Lw/(N1-1)
h2    = Lh/(N2-1)
nunk  = N1*N2

gamma = (2*mu)*(1/(h1**2) + 1/(h2**2))
alpha = - mu/(h1**2)
beta  = - mu/(h2**2)

def AssemblyMatDiags(N1,N2,nunk,alpha,beta,gamma):
   d1 = np.ones(nunk)    * gamma
   d2 = np.ones(nunk-1)  * alpha
   d3 = np.ones(nunk-N1) * beta

   d1[range(0,N1)]                = 1
   d1[range(N1*(N2-1),nunk)]      = 1
   d1[range(N1,N1*(N2-1)+1,N1)]   = 1

   d2[range(0,N1)]                = 0
   d2[range(N1*(N2-1),nunk-1)]    = 0
   d2[range(N1,N1*(N2-1)+1,N1)]   = 0
   d2[range(2*N1-1,N1*(N2-1),N1)] = 0

   d3[range(0,N1)]                 = 0
   d3[range(N1*(N2-1),nunk-N1)]    = 0
   d3[range(N1,N1*(N2-1)+1-N1,N1)] = 0
   d3[range(2*N1-1,N1*(N2-1),N1)]  = 0

   d3i = d3[::-1]
   d2i = d2[::-1]
   return diags([d3i,d2i,d1,d2,d3], [-N1,-1,0,1,N1])

def AssemblyMatb(N1,N2,h1,h2,nunk):
   b                              = np.ones((nunk,1))
   b[range(0,N1)]                 = 0
   b[range(N1*(N2-1),N1*N2)]      = 0
   b[range(N1,N1*(N2-1)+1,N1)]    = 0
   b[range(2*N1-1,N1*(N2-1),N1)]  = 0
   bm                             = np.ones((nunk,1)) * h1 * h2
   bm[range(0,N1)]                = 0
   bm[range(N1*(N2-1),N1*N2)]     = 0
   bm[range(N1,N1*(N2-1)+1,N1)]   = 0
   bm[range(2*N1-1,N1*(N2-1),N1)] = 0
   return b , bm

### Linear system
A      = AssemblyMatDiags(N1,N2,nunk,alpha,beta,gamma)
b, bm  = AssemblyMatb(N1,N2,h1,h2,nunk)
Id     = np.identity(nunk)
Aleft  = sparse.csr_matrix((ro/dt)*Id +    theta *A)
Aright = sparse.csr_matrix((ro/dt)*Id - (1-theta)*A)

### Declaration
u  = np.zeros((nunk,1))
un = np.zeros((nunk,1))
q  = np.zeros((nt+1,1))
gg = np.zeros((nt+1,1))

q[0]  = np.dot(np.transpose(bm),un)
gg[0] = (1-theta)*np.cos(2*np.pi*tn/T) + theta*np.cos(2*np.pi*(tn+dt)/T)

### Creating coordinates
x = np.linspace(0.0, Lw, N1)
y = np.linspace(0.0, Lh, N2)
z = np.linspace(0.0, Lz, N3)

def savevtkfile(namearquive, u, x, y, z, N1, N2, N3):
    with open(namearquive, "w") as f:
       header =  "# vtk DataFile Version 3.0\n"
       header += "vtk output\n"
       header += "ASCII\n"
       header += "DATASET STRUCTURED_GRID\n"
       header += "DIMENSIONS "+str(N1)+" "+str(N2)+" "+str(N3)+"\n"
       header += "POINTS "+str(N1*N2*N3)+" FLOAT\n"
       f.write(header)
       for k in range(N3):
          for j in range(N2):
             for i in range(N1):
                f.write('%.16f %.16f %.16f\n' % (x[i],y[j],z[k]))
       f.write("POINT_DATA "+str(N1*N2*N3)+" \n")
       #f.write("\n")
       f.write("VECTORS u FLOAT\n")
       for k in range(N3):
          for j in range(N2):
             for i in range(N1):
                f.write('%.18f %.18f %.18f\n' % (0.0, 0.0, U[i][j]))

####### Loop TIME STEPS ########
i   = 1
while (i <= nt):
   ### Source term
   gg[i] = (1-theta)*np.cos(2*np.pi*tn/T) + theta*np.cos(2*np.pi*(tn+dt)/T)
   ### Sistema [A*un = gg]
   u     = spsolve(Aleft, Aright*un - gg[i]*b)
   ### Reshape u for solution
   u     = u.reshape((nunk,1))
   ### Copy u to un
   un    = np.copy(u)
   
   q[i]  = np.dot(np.transpose(bm), u)
   tn    += dt
   
   if (i % 100 == 0):
      print('Iteracao := ',i,' de ',nt)
      namearquive = 'grid'+str(i).zfill(len(str(nt)))+'.vtk'
      U     = np.reshape(u, (N1, N2))
      ### 2D
      savevtkfile(namearquive, U, x, y, z, N1, N2, 1)
      ### 3D
      #savevtkfile(namearquive, U, x, y, z, N1, N2, N3)
      shutil.move(cwd+'/'+namearquive,dir_)
   
   i   += 1
 
np.savetxt(dir_+'/u.txt', u)
np.savetxt(dir_+'/q.txt', q)
np.savetxt(dir_+'/g.txt', gg)
    
# tt = np.linspace(0, tf, nt)
# plt.plot(tt,q/(-max(q)), tt,gg)
# plt.legend(['$\\frac{Q}{-max(Q)}$', 'G(t)'])
# plt.xlabel('t')
# plt.title('$\\mu=$',mu)
# plt.savefig()
