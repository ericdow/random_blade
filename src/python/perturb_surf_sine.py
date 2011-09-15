import read_blade
import simulate
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

rpath = '/home/ericdow/code/random_blade/input/blade_surf.dat'
wpath = '/home/ericdow/code/random_blade/input/blade_surf_mod.dat'

cdim, sdim, x, y, z = read_blade.read_coords(rpath)
xu, yu, zu, xl, yl, zl = read_blade.split_blade(cdim, sdim, x, y, z)
xmc = 0.5*(xu+xl)
ymc = 0.5*(yu+yl)
zmc = 0.5*(zu+zl)

dx = xu - xmc
dy = yu - ymc
dz = zu - zmc

# camber random process
sc,tc = read_blade.xyz2st(xmc,ymc,zmc)
ns = xmc.shape[0]
nt = xmc.shape[1]
Ks = int(ns/5.)          # mode cutoff
Kt = int(nt/5.)
M = 1
N = 1
fc = 0.0005*sc[-1]*outer(sin(2*pi*M*sc/sc[-1]),sin(2*pi*N*tc/tc[-1]))

tg, sg = meshgrid(tc,sc)
                     
nc = read_blade.calcNormals(xmc,ymc,zmc)

xmc = xmc + nc[:,:,0]*fc
ymc = ymc + nc[:,:,1]*fc
zmc = zmc + nc[:,:,2]*fc
 
xu = xmc + dx
yu = ymc + dy
zu = zmc + dz

xl = xmc - dx
yl = ymc - dy
zl = zmc - dz

xp = vstack((xu[:-1,:],xl[::-1,:]))
yp = vstack((yu[:-1,:],yl[::-1,:]))
zp = vstack((zu[:-1,:],zl[::-1,:]))

# normal random process
sp,tp = read_blade.xyz2st(xp,yp,zp)
ns = xp.shape[0]
nt = xp.shape[1]
Ks = int(ns/4.)
Kt = int(nt/4.)

# fp = simulate.randProcessPeriodic(sp, tp, ns, nt, Ks, Kt)
M = 2
N = 1
fp = 0.05*sp[-1]*outer(sin(2*pi*M*sp/sp[-1]),cos(2*pi*N*tp/tp[-1]))

######
# MAY BE EASIER TO JUST INTERPOLATE THE K-L EIGENFUNCTIONS
# IN 1D AND RECONSTRUCT FROM THAT
######
'''
# interpolate back to original mesh
ss = linspace(0,sp[-1]-sp[0],ns)
tt = linspace(0,tp[-1]-tp[0],nt)
ff = simulate.randProcessPeriodic(ss, tt, ns, nt, Ks, Kt)
fp = simulate.interp(ss, tt, ff, sp, tp)

tg, sg = meshgrid(tt,ss)
pylab.contourf(sg,tg,ff,30)
pylab.contourf(sg+sg[-1],tg,ff,30)
pylab.colorbar()
pylab.axes().set_aspect('equal', 'datalim')
pylab.figure()
tg, sg = meshgrid(tp,sp)
pylab.contourf(sg,tg,fp,30)
pylab.contourf(sg+sg[-1],tg,fp,30)
pylab.colorbar()
pylab.axes().set_aspect('equal', 'datalim')
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

'''
fig = pylab.figure()
ax = Axes3D(fig)
# surf = ax.plot_surface(xp, yp, zp, rstride=1, cstride=1, cmap = cm.jet)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.jet(fp),\
                       linewidth=0, antialised=0, shade=False)
pylab.show()
'''


