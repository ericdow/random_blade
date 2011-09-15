import read_blade
import simulate
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

path = '../../input/blade_surf_mod_proj.dat'

lines = file(path).readlines()

fdim = 30
sdim = 113
 
x = array([line.strip().split()[0] for line in lines], float).T
y = array([line.strip().split()[1] for line in lines], float).T
z = array([line.strip().split()[2] for line in lines], float).T
nc = x.size

x = reshape(x, (sdim,fdim))
y = reshape(y, (sdim,fdim))
z = reshape(z, (sdim,fdim))
 
fig = pylab.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap = cm.jet)
pylab.show()
