import read_blade
import simulate
from math import *
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
##########
# SINE
fc = 0.0010*sc[-1]*outer(sin(2*pi*M*sc/sc[-1]),sin(2*pi*N*tc/tc[-1]))
##########
'''
##########
# ANGLE
zmid = 0.5*(max(zmc[:,0]) + min(zmc[:,0]))
chord = sqrt((max(zmc[:,0])-min(zmc[:,0]))**2 + (max(zmc[:,0])-min(zmc[:,0]))**2)
# fc = 0.005*chord*(2.0*(zmc-min(zmc[:,0]))/(max(zmc[:,0])-min(zmc[:,0]))-1.0)
fc = 0.01*chord*(1.0 - zmc/max(zmc[:,0]))
dx1 = xmc[0,0]-xmc[-1,0]
dy1 = ymc[0,0]-ymc[-1,0]
dz1 = zmc[0,0]-zmc[-1,0]
##########
'''
'''
##########
# TWIST WITH TIP FIXED
chord = sqrt((max(zmc[:,0])-min(zmc[:,0]))**2 + (max(zmc[:,0])-min(zmc[:,0]))**2)
fc = -0.01*chord*(1.0 - zmc/max(zmc[:,0]))
for i in range(ns):
    rh = sqrt(xmc[i,0]**2 + ymc[i,0]**2)
    rt = sqrt(xmc[i,-1]**2 + ymc[i,-1]**2) - 0.00356
    for j in range(nt):
        r = sqrt(xmc[i,j]**2 + ymc[i,j]**2)
        fc[i,j] *= 1.0 - (r-rh)/(rt-rh)
dx1 = xmc[0,0]-xmc[-1,0]
dy1 = ymc[0,0]-ymc[-1,0]
dz1 = zmc[0,0]-zmc[-1,0]
##########
'''

tg, sg = meshgrid(tc,sc)
                     
nc = read_blade.calcNormalsCamber(xmc,ymc,zmc)

xmc = xmc + nc[:,:,0]*fc
ymc = ymc + nc[:,:,1]*fc
zmc = zmc + nc[:,:,2]*fc

# pl = -1
# pylab.figure()
# pylab.plot(zmc[:,pl],ymc[:,pl],'b*-')
# pylab.plot(z[:,pl],y[:,pl])
# pylab.plot(zmc[:,pl],ymc[:,pl],'r*-')
 
dx2 = xmc[0,0]-xmc[-1,0]
dy2 = ymc[0,0]-ymc[-1,0]
dz2 = zmc[0,0]-zmc[-1,0]

# dth = acos(1./sqrt(dx1**2+dy1**2+dz1**2)/sqrt(dx2**2+dy2**2+dz2**2)*(dx1*dx2+dy1*dy2+dz1*dz2))*180/pi
# print 'ANGLE CHANGE: ', dth

xu = xmc + dx
yu = ymc + dy
zu = zmc + dz

xl = xmc - dx
yl = ymc - dy
zl = zmc - dz

xp = vstack((xu[:-1,:],xl[::-1,:]))
yp = vstack((yu[:-1,:],yl[::-1,:]))
zp = vstack((zu[:-1,:],zl[::-1,:]))

# pylab.plot(zp[:,pl],yp[:,pl])
# pylab.axis('equal')

# write out the blade surface
f = open(wpath,'w')
f.write('CDIM:       %d\n' % cdim)
f.write('SDIM:       %d\n' % sdim)
for i in arange(cdim):
    for j in arange(sdim):
        f.write('%20.8E' * 3 % (xp[i,j],yp[i,j],zp[i,j]))
        f.write('\n')

f.close()

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
# surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.jet(fp),\
#                        linewidth=0, antialised=0, shade=False)
# pylab.show()
