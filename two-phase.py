import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# permeability
k1 = 1.0e-3
k2 = 1.0e-5
# viscosity
mu_w, mu_o = 1.0, 5.0
# porosity
phi = 0.2

# pressure boundary conditions
aw = 1.0
pw_inj = 1.2
pw_prod = 1.0

# saturation initial and boundary condition
sinit = 0.2
sw_inj = 1.0
# islin = True
islin = False

def lmbw(s):
    if islin: return s/mu_w
    return s*s/mu_w
def lmbn(s):
    if islin: return (1-s)/mu_o
    return (1-s)*(1-s)/mu_o
def lmb(s):
    return (lmbn(s) + lmbw(s))
def fw(s):
    return lmbw(s)/(lmbn(s) + lmbw(s))

def getIndex(i, j, Nx):
    return j*Nx+i

def plot_solution(p,s, filename):
    fig = plt.figure(figsize=(9, 3))
    ax1 = fig.add_subplot(121)
    im1 = ax1.imshow(np.reshape(p, (Ny, Nx)),  cmap=plt.get_cmap('jet'))
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='5%', pad=0.1)
    fig.colorbar(im1, cax=cax1, orientation='vertical')
    ax2 = fig.add_subplot(122)
    im2 = ax2.imshow(np.reshape(s, (Ny, Nx)),  cmap=plt.get_cmap('jet'))
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes('right', size='5%', pad=0.1)
    fig.colorbar(im2, cax=cax2, orientation='vertical')
#     plt.show()
    plt.savefig(filename)
    plt.close('all')

# ----- generate permeability, k(x) -----
Nx, Ny = 20, 20
arrK = np.zeros(Nx*Ny)
for i in range(Nx):
    for j in range(Ny):
        I = getIndex(i, j, Nx)
        if j > Ny*2/3 and i < Nx/2:
            arrK[I] = k1
        else:
            arrK[I] = k2
plot_solution(arrK, arrK, 'k.png')


# time parameters
Tcount = 100
Tmax = 1.0e5
Tsave = int(Tcount/5)
dt = Tmax/Tcount

hx, hy = 1.0/Nx, 1.0/Ny
volK = hx*hy

A = np.zeros((Nx*Ny, Nx*Ny))
b = np.zeros(Nx*Ny)
r = np.zeros(Nx*Ny)

snew = np.zeros(Nx*Ny)
s = np.zeros(Nx*Ny)
p = np.zeros(Nx*Ny)

# initial condition
for i in range(Nx):
    for j in range(Ny):
        I = getIndex(i, j, Nx)
        s[I] = sinit

t = 0.0
tcounter = 0
while (tcounter < Tcount+1):
    # ----- pressure system -----
    for i in range(Nx):
        for j in range(Ny):
            I = getIndex(i, j, Nx)
            diagval = 0
            if j!=0:
                J = getIndex(i, j-1, Nx)
                sval = (s[I] + s[J])/2 # average
                kval = 2.0/(1.0/arrK[J]+1.0/arrK[I]) # harmonic average
                val = lmb(sval)*kval*(hx/hy)
                A[I,J] = -val; diagval += val
            if j!=(Ny-1):
                J = getIndex(i, j+1, Nx)
                sval = (s[I] + s[J])/2 # average
                kval = 2.0/(1.0/arrK[J]+1.0/arrK[I]) # harmonic average
                val = lmb(sval)*kval*(hx/hy)
                A[I,J] = -val; diagval += val
            if i!=0:
                J = getIndex(i-1, j, Nx)
                sval = (s[I] + s[J])/2 # average
                kval = 2.0/(1.0/arrK[J]+1.0/arrK[I]) # harmonic average
                val = lmb(sval)*kval*(hy/hx)
                A[I,J] = -val; diagval += val
            if i!=(Nx-1):
                J = getIndex(i+1, j, Nx)
                sval = (s[I] + s[J])/2 # average
                kval = 2.0/(1.0/arrK[J]+1.0/arrK[I]) # harmonic average
                val = lmb(sval)*kval*(hy/hx)
                A[I,J] = -val; diagval += val
            # Robyn boundary conditons
            if i == 0 or i == (Nx-1):
                if i == 0:
                    sval = sw_inj
                    pval = pw_inj
                if i == (Nx-1):
                    sval = s[I]
                    pval = pw_prod
                kval = arrK[I]
                val = aw*lmb(sval)*kval*(hy/hx/2)
                b[I] = val*pval
                r[I] = -val; diagval += val
            # set diagonal
            A[I,I] = diagval
    # ----- pressure solve -----
    p = np.linalg.solve(A, b)
    
    # ----- saturation explicit -----
    for i in range(Nx):
        for j in range(Ny):
            I = getIndex(i, j, Nx)
            # all flows
            sums = 0.0
            if j!=0:
                J = getIndex(i, j-1, Nx)
                vel = -A[I,J]*(p[I] - p[J])
                if vel < 0:  
                    sval = s[J]
                else:
                    sval = s[I]
                sums += fw(sval)*vel  
            if j!=(Ny-1):
                J = getIndex(i, j+1, Nx)
                vel = -A[I,J]*(p[I] - p[J])
                if vel < 0:  
                    sval = s[J]
                else:
                    sval = s[I]
                sums += fw(sval)*vel 
            if i!=0:
                J = getIndex(i-1, j, Nx)
                vel = -A[I,J]*(p[I] - p[J])
                if vel < 0:  
                    sval = s[J]
                else:
                    sval = s[I]
                sums += fw(sval)*vel 
            if i!=(Nx-1):
                J = getIndex(i+1, j, Nx)
                vel = -A[I,J]*(p[I] - p[J])
                if vel < 0:  
                    sval = s[J]
                else:
                    sval = s[I]
                sums += fw(sval)*vel 
            # boundary condition
            if i == 0 or i == (Nx-1):
                if i == 0:
                    sval = sw_inj
                    pval = pw_inj
                if i == (Nx-1):
                    sval = s[I]
                    pval = pw_prod
                vel = -r[I]*(p[I] - pval) 
                sums += fw(sval)*vel
            # explicit calculation of saturation
            ss = s[I] - dt/(volK*phi)*sums
            ss = max(min(ss, sw_inj), sinit)
            snew[I] = ss
    # update
    for I in range(Nx*Ny):
        s[I] = snew[I]
    
    t += dt
    tcounter += 1
    
    if tcounter%Tsave == 0:
        plot_solution(p, s, 'sol-t'+str(tcounter)+'.png')
        print('Time[%d] = %g, s in [%g, %g]' % (tcounter, t, s.min(), s.max()))
    