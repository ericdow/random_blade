import pylab, os, math, string
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import read_blade

inp = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/'
rundir = '/mnt/pwfiles/ericdow/random_blade/runs/'

# read in original blade surface
cdim, sdim, xn, yn, zn = read_blade.read_coords(inp+'blade_surf.dat')

# read in blade surfaces of completed runs
iran = zeros((len(os.listdir(rundir))))
i = 0
for rd in os.listdir(rundir):
    if os.path.exists(rundir+str(rd)+'/utcfd_out.monitoring.dat'):
        lines = file(rundir+str(rd)+'/utcfd_out.monitoring.dat').readlines()
        if lines[-1][1:4] == 'Job':
            iran[i] = 1
        lines = file(rundir+str(rd)+'/Z.dat').readlines()
        npca = len(lines)
    i += 1

nruns = int(sum(iran))

# perturbed coordinates
xp = zeros((nruns,cdim,sdim))
yp = zeros((nruns,cdim,sdim))
zp = zeros((nruns,cdim,sdim))

ir = 0
i = 0
for rd in os.listdir(rundir):
    if iran[ir]:
        c,s,xp[i,:,:],yp[i,:,:],zp[i,:,:] = read_blade.read_coords(rundir+str(rd)+'/blade_surf_mod.dat')
        i += 1
    ir += 1

# plot a few blade sections
pylab.plot(xn[:,0],yn[:,0],linewidth=2)
for i in range(10):
    pylab.plot(xp[i,:,0],yp[i,:,0])
pylab.title('Hub')
pylab.axis('Equal')

pylab.figure()
pylab.plot(xn[:,sdim/2],yn[:,sdim/2],linewidth=2)
for i in range(10):
    pylab.plot(xp[i,:,sdim/2],yp[i,:,sdim/2])
pylab.title('Midspan')
pylab.axis('Equal')

pylab.figure()
pylab.plot(xn[:,sdim-1],yn[:,sdim-1],linewidth=2)
for i in range(10):
    pylab.plot(xp[i,:,sdim-1],yp[i,:,sdim-1])
pylab.title('Shroud')
pylab.axis('Equal')

pylab.show()
