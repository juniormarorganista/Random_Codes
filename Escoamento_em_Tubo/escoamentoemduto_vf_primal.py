import numpy as np
from numpy import linalg
import scipy.sparse
import scipy.sparse.linalg
import time
import os
import shutil
import matplotlib
import matplotlib.pyplot as plt
from uvw import RectilinearGrid, DataArray

from scipy.sparse.linalg import inv

# IO
cwd = os.getcwd()
dir_=cwd+'/results_U'

if os.path.isdir(dir_):
    # remove dir_ and all contains
    shutil.rmtree(dir_)

if not os.path.exists(dir_):
    os.makedirs(dir_)

print(dir_)

#True = Activate this for efficiency on large grids (Com fluxo e mais rapido)
flag_coo = False 
#OBS  = Para calcular o sistema primal somente com U, use False nas duas flags acima
#True = Mu constante em todo dominio e False = mu interpolate(item 5)
flag_mu  = False  

# Parameters
mu0 = 1.e-3
rho = 1.e3

# Creating coordinates
L1, N1 = 20.e-6, 101
L2, N2 = 10.e-6, 101
n1, n2 = N1-1, N2-1
Nc     = n1*n2
Nh     = n1 * (n2+1)
Nv     = n2 * (n1+1)
Nf     = Nh + Nv

nunk = Nc 

def CreateMesh(L1, L2, N1, N2, refine_walls=False, alpha=0.12):
    x = np.linspace(0.0, L1, N1)
    y = np.linspace(0.0, L2, N2)
    if(False):
        for i1 in range(1,N1-1):
            x[i1] = alpha*(x[i1]-0.5)
            x[i1] = 0.5 * (np.cos(np.pi * (x[i1] - 1.) / 2.) + 1.)
        
        for i2 in range(1,N2-1):
            y[i2] = alpha*(y[i2]-0.5)
            y[i2] = 0.5 * (np.cos(np.pi * (y[i2] - 1.) / 2.) + 1.)

    if(refine_walls):
        for i1 in range(1,N1-1):
            x[i1] = x[i1] - alpha*np.sin(2*np.pi*x[i1])
        
        for i2 in range(1,N2-1):
            y[i2] = y[i2] - alpha*np.sin(2*np.pi*y[i2])

    return x, y

x, y = CreateMesh(L1, L2, N1, N2, False, 0.12) # Be careful with alpha

# Post-processing
def WriteSol(x, y, n1, n2, nstep, n, tn, U):
    filename = 'ref3FVsol_primal_'+str(n).zfill(len(str(nstep)))+'.vtr'
    grid = RectilinearGrid(filename, (x, y), compression=True)
    grid.addCellData(DataArray(U.reshape(n1,n2), range(2), 'Velocity'))
    grid.write()
    shutil.move(cwd+'/'+filename ,dir_)

# Just for testing the mesh
#U = np.random.rand(n1*n2)
#WriteSol(x, y, n1, n2, 0, 0.0, U)

def jglob(i1,i2,n1):
    return i1 + i2*n1

def fh(i1,i2,n1):
    return i1 + i2*n1

def fv(i1,i2,n1,n2):
    return n1*(n2+1) + i1 + i2*(n1+1)

def Visc(x, y, L1, L2, mu0):
    pot = -50*( ((x-0.5*L1)**2)/L1**2 + ((y-0.5*L2)**2)/L2**2 )
    mu  = mu0 + mu0*np.exp(pot)
    return mu

# Face-to-cell connectivity
FC = -np.ones(shape=(Nf,2), dtype=int)
for k1 in range(n1):
    for k2 in range(n2+1):
        f = fh(k1,k2,n1)
        if(k2 == 0):
            FC[f,0] = jglob(k1,k2,n1)
            FC[f,1] = -1
        elif(k2 == n2):
            FC[f,0] = jglob(k1,k2-1,n1)
            FC[f,1] = -3
        else:
            FC[f,0] = jglob(k1,k2-1,n1)
            FC[f,1] = jglob(k1,k2,n1)

for k1 in range(n1+1):
    for k2 in range(n2):
        f = fv(k1,k2,n1,n2)
        if(k1 == 0):
            FC[f,0] = jglob(k1,k2,n1)
            FC[f,1] = -4
        elif(k1 == n1):
            FC[f,0] = jglob(k1-1,k2,n1)
            FC[f,1] = -2
        else:
            FC[f,0] = jglob(k1-1,k2,n1)
            FC[f,1] = jglob(k1,k2,n1)

# Cell-to-face connectivity
CF = np.zeros(shape=(Nc,4), dtype=int)
for i1 in range(n1):
    for i2 in range(n2):
        g       = jglob(i1,i2,n1)
        CF[g,:] = [fh(i1,i2,n1), fv(i1+1,i2,n1,n2), fh(i1,i2+1,n1), fv(i1,i2,n1,n2)]

#plt.spy(A.todense(),precision=0.1,markersize=1)
def BuildSystem(nunk, x, y, CF, FC, mu0, rho, theta, Deltat):
    #------------------------------------------------
    # Total number of unknowns
    #------------------------------------------------
    N1, N2 = len(x), len(y)
    n1, n2 = N1-1, N2-1
    Nc = n1 * n2
    Nh = n1 * (n2+1)
    Nv = n2 * (n1+1)
    Nf = Nh + Nv
    if(flag_coo):
        row, col, coefL, coefR = [], [], [], []
    #------------------------------------------------
    # Preliminaries
    #------------------------------------------------
    t0 = time.time()
    vcells = np.zeros(Nc)
    Dx = np.zeros(Nc)
    Dy = np.zeros(Nc)
    xc = np.zeros(Nc)
    yc = np.zeros(Nc)
    for i1 in range(n1):
        for i2 in range(n2):
            g = jglob(i1,i2,n1)
            Dx[g] = (x[i1+1] - x[i1])
            Dy[g] = (y[i2+1] - y[i2])
            xc[g] = x[i1] + 0.5*Dx[g]
            yc[g] = y[i2] + 0.5*Dy[g]
    vcells = Dx * Dy
    #------------------------------------------------
    # Matrix K and Kinv
    #------------------------------------------------
    dgnKinv = np.zeros(Nf)
    dgnK    = np.zeros(Nf)
    fs      = np.zeros(Nf)
    t0      = time.time()
    for f in range(Nf):
        g, n = FC[f, :]
        ### Calculate mu
        if (flag_mu):
            mu = mu0
        else:
            if(n*g > 0):
                if (n-g == 1):
                    mu = Visc(xc[g] + 0.5*Dx[g], yc[g], L1, L2, mu0)
                else:
                    mu = Visc(xc[g], yc[g] + 0.5*Dy[g], L1, L2, mu0)
            else:
                if (n == -1):
                    mu = Visc(xc[g], yc[g] - 0.5*Dy[g], L1, L2, mu0)
                elif(n == -2):
                    mu = Visc(xc[g] + 0.5*Dx[g], yc[g], L1, L2, mu0)
                elif(n == -3):
                    mu = Visc(xc[g], yc[g] + 0.5*Dy[g], L1, L2, mu0)
                elif(n == -4):
                    mu = Visc(xc[g] - 0.5*Dx[g], yc[g], L1, L2, mu0)
        ### Diagonals K and Kinv
        fs[f] = Dx[g] if (f < Nh) else Dy[g]
        if(g >= 0 and n >= 0):
            dv         = np.array([xc[g]-xc[n], yc[g]-yc[n]])
            dgnKinv[f] = np.linalg.norm(dv) / mu
            dgnK[f]    = mu / np.linalg.norm(dv)
        else:
            dgnKinv[f] = (0.5*Dy[g]) / mu if(f < Nh) else (0.5*Dx[g]) / mu
            dgnK[f]    = mu / (0.5*Dy[g]) if(f < Nh) else mu / (0.5*Dx[g])
        if(flag_coo):
            row.append(f)
            col.append(f)
            coefL.append(dgn[f]*theta)
            coefR.append(dgn[f]*(theta-1.0))

    Kinv = scipy.sparse.diags([dgnKinv], [0], format='lil')
    K    = scipy.sparse.diags([dgnK], [0], format='lil')
    print('Matrix K: ', time.time() - t0)
    #------------------------------------------------
    # Matrix A
    #------------------------------------------------
    t0 = time.time()
    A = scipy.sparse.lil_matrix((Nf,Nc), dtype=float)
    for f in range(Nf):
        g, n = FC[f, :]
        if(g >= 0 and n >= 0):
            A[f,g] = +1.0
            A[f,n] = -1.0
            if(flag_coo):
                row.append(f)
                row.append(f)
                col.append(Nf+g)
                col.append(Nf+n)
                coefL.append(-1.0*theta)
                coefL.append(1.0*theta)
                coefR.append(-1.0*(theta-1.0))
                coefR.append(1.0*(theta-1.0))
        else:
            A[f,g] = +1.0
            if(flag_coo):
                row.append(f)
                col.append(Nf+g)
                coefL.append(-1.0*theta)
                coefR.append(-1.0*(theta-1.0))
    print('Matrix A: ', time.time() - t0)
    #------------------------------------------------
    # Matrix C
    #------------------------------------------------
    t0 = time.time()
    C  = scipy.sparse.lil_matrix((Nc,Nf), dtype=float)
    for g in range(Nc):
        for f in CF[g,:]:
            C[g,f] = fs[f] if(g == FC[f,0]) else -fs[f]
            if(flag_coo):
                row.append(Nf+g)
                col.append(f)
                coefL.append(C[g,f]*theta)
                coefR.append(C[g,f]*(theta-1.0))
    print('Matrix C: ', time.time() - t0)
    #------------------------------------------------
    # Matrix M
    #------------------------------------------------
    t0 = time.time()
    d = (rho/Deltat) * vcells
    M = scipy.sparse.diags([d], [0], format='lil')
    if(flag_coo):
        for g in range(Nc):
            row.append(Nf+g)
            col.append(Nf+g)
            coefL.append(M[g,g])
            coefR.append(M[g,g])
    print('Matrix M: ', time.time() - t0)
    #------------------------------------------------
    # Global matrix
    #------------------------------------------------
    if(flag_coo):
        t0    = time.time()
        row   = np.array(row)
        col   = np.array(col)
        coefL = np.array(coefL)
        coefR = np.array(coefR)
        MatL  = scipy.sparse.coo_matrix((coefL, (row, col)), shape=(nunk, nunk))
        MatR  = scipy.sparse.coo_matrix((coefR, (row, col)), shape=(nunk, nunk))
        print('Matrix L and R (coo format): ', time.time() - t0)
        # Vector gb
        gb = np.zeros(nunk)
        gb[Nf:nunk] = vcells[0:Nc]
    else:
        t0      = time.time()
        CKA     = C @ K @ A
        gb      = np.zeros(Nc)
        gntheta = np.zeros(Nf)
        MatL    = M + theta * CKA 
        MatR    = (theta - 1) * ( CKA + M )
        gb      = np.zeros(nunk)
        gb      = vcells[0:Nc]  - C @ K @ gntheta
        print('Matrix L and R (csr-slicing) - only U: ', time.time() - t0)
    #------------------------------------------------
    #Visualize sparsity pattern
    #plt.spy(B.todense(),precision=0.1,markersize=1)
    #plt.show()
    #------------------------------------------------
    return MatL, MatR, gb

def PressGrad(t):
    T = 2.0
    return np.cos(2.0*np.pi*t/T)

#Time interval
t0, tf   = 0.0, 1.0
nsteps   = 1000
theta    = 0.5
freq_out = 100

#Initialization
t, Deltat = np.linspace(t0, tf, nsteps, retstep=True, endpoint=False)

# Build algebraic objects
MatL, MatR, gb = BuildSystem(nunk, x, y, CF, FC, mu0, rho, theta, Deltat)

# Initial condition
U = np.zeros(nunk)

print('Begin loop over time steps:')
t0 = time.time()
for n in range(nsteps):
    if(n % freq_out == 0 or n == nsteps-2):
        print('Step solution:= ', n, 'at time := ', t[n])
    #    print('Writing solution file', n, 'at time= ', t[n])
    #    if (flag_JU):
    #        U = JU[Nf:nunk]
    #    else:
    #        U=JU
    #    WriteSol(x, y, n1, n2, n, t[n], U)
    G   = PressGrad(t[n] + theta*Deltat)
    rhs = MatR @ U - G * gb
    U  = scipy.sparse.linalg.spsolve(MatL, rhs)

    WriteSol(x, y, n1, n2, nsteps, n, t[n], U)

print('Time stepping: ', time.time()-t0)
