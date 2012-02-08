import read_blade
import simulate
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

rpath = '/home/ericdow/code/random_blade/input/blade_surf.dat'
wpath = '/home/ericdow/code/random_blade/input/blade_surf_mod.dat'

cdim, sdim, x, y, z = read_blade.read_coords(rpath)

xp = x
yp = y
zp = z

# normal random process
sp,tp = read_blade.xyz2st(xp,yp,zp)
ns = xp.shape[0]
nt = xp.shape[1]
Ks = ns
Kt = nt
M = 1
N = 1
fp = 0.0010*sp[-1]*outer(cos(2*pi*M*sp/sp[-1]),cos(2*pi*N*tp/tp[-1]))
fp[0,:] = fp[-1,:]

'''
tg, sg = meshgrid(tp,sp)
pylab.contourf(sg,tg,fp,30)
# pylab.contourf(sg+sg[-1],tg,fp,30)
pylab.colorbar()
pylab.xlim(0.0,sp[-1])
pylab.axes().set_aspect('equal', 'datalim')
pylab.xlabel('x')
pylab.ylabel('y')
pylab.show()
'''

np = read_blade.calcNormals(xp,yp,zp)

xp = xp + np[:,:,0]*fp
yp = yp + np[:,:,1]*fp
zp = zp + np[:,:,2]*fp

# write out the blade surface
f = open(wpath,'w')
f.write('CDIM:       %d\n' % cdim)
f.write('SDIM:       %d\n' % sdim)
for i in arange(cdim):
    for j in arange(sdim):
        f.write('%20.8E' * 3 % (xp[i,j],yp[i,j],zp[i,j]))
        f.write('\n')

f.close()

# write out the normal perturbation
g = open('normal_field.dat','w')
for i in arange(cdim):
    for j in arange(sdim):
        g.write('%20.8E' * 4 % (x[i,j],y[i,j],z[i,j],fp[i,j]))
        g.write('\n')

g.close()

'''
fig = pylab.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(sg, tg, fp, rstride=1, cstride=1, cmap = cm.jet)
'''

'''
fig = pylab.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(tg, sg, nf[1:-1,1:-1,2], rstride=1, cstride=1, cmap = cm.jet)
# fig.colorbar(surf)
'''

fig = pylab.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(xp, yp, zp, rstride=1, cstride=1, cmap = cm.jet)
#surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.jet(fp),\
#                       linewidth=0, antialised=0, shade=False)
pylab.show()
       

