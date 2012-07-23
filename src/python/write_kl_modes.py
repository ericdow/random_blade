import read_blade, simulate, write_tecplot, project_GP
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# where to read in blade geometry
rpath = 'blade_surf.dat'

# number of modes to include
nkl = 50

cdim, sdim, x, y, z = read_blade.read_coords(rpath)

# generate K-L modes
s,t = read_blade.xyz2st(x,y,z)
ns = x.shape[0]
nt = x.shape[1]
Ks = ns
Kt = nt
f, SX, UX, SY, UY = simulate.randProcessPeriodic(s, t, ns, nt, Ks, Kt)

# sort the K-L modes in order of eigenvalue
eigs = outer(SX,SY)
eigs_1d = reshape(eigs,ns*nt,1)
ind = flipud(argsort(eigs_1d)[-nkl:])
ind_2d = zeros((nkl,2))
for i in range(nkl):
    ind_2d[i,:] = unravel_index(ind[i],(nt,ns))
ind_2d = fliplr(ind_2d)

# write out the K-L modes/eigenvalues
f = open('kl_modes/S.dat','w')
for n in range(nkl):
    f.write('%e\n' % eigs_1d[ind[n]])
f.close()
for n in range(nkl):
    imode = ind_2d[n,0]
    jmode = ind_2d[n,1]
    V = outer(UX[:,imode],UY[:,jmode])
    # project K-L mode onto mesh geometry
    sABCD = array([0.05, 0.45, 0.55, 0.95])
    slope = 0.1 # how much stuff is pushed towards hub/shroud
    V_proj = project_GP.project(s,t,sABCD,V,x,y,z,slope)

    f = open('kl_modes/V'+str(n+1)+'.dat','w')
    f.write('CDIM:       %d\n' % cdim)
    f.write('SDIM:       %d\n' % sdim)
    for i in arange(cdim):
        for j in arange(sdim):
            f.write('%e\n' % V_proj[i,j])
    f.close()

    write_tecplot.write_blade_surf(x,y,z,V_proj,'blade.dat') 
