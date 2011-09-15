import read_blade
import simulate
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time

rpath = '/home/ericdow/code/random_blade/input/blade_surf.dat'
wpath = '/home/ericdow/code/random_blade/input/blade_surf_mod.dat'

cdim, sdim, x, y, z = read_blade.read_coords(rpath)
xu, yu, zu, xl, yl, zl = read_blade.split_blade(cdim, sdim, x, y, z)

xp = vstack((xu[:-1,:],xl[::-1,:]))
yp = vstack((yu[:-1,:],yl[::-1,:]))
zp = vstack((zu[:-1,:],zl[::-1,:]))

sp,tp = read_blade.xyz2st(xp,yp,zp)
ns = xp.shape[0]
nt = xp.shape[1]
Ks = int(ns/5.)
Kt = int(nt/5.)

# 1D test
i1 = 10
i2 = 50

tp = linspace(0,tp[-1],nt)

N = 20000
f_sample = zeros((2,N))
for n in arange(N):
    tic = time.clock()
    f = simulate.randProcess1DPeriodic(tp, nt, nt)
    toc = time.clock()
    print n,':',toc-tic
    f_sample[0,n] = f[i1]
    f_sample[1,n] = f[i2]

mu1 = sum(f_sample[0,:])/N
mu2 = sum(f_sample[1,:])/N

cov = dot(f_sample[0,:]-mu1,f_sample[1,:]-mu2)/N

print 'estimate:',cov

th1 = tp[i1]/tp[-1]*2*pi
th2 = tp[i2]/tp[-1]*2*pi
d = sqrt((cos(th1)-cos(th2))**2 + (sin(th1)-sin(th2))**2)
cov_ex = simulate.covfun(d)
print 'exact:', cov_ex

mu = zeros(N)
for n in arange(N):
    mu[n] = sum(f_sample[0,:n+1])
    mu[n] /= (n+1)
pylab.loglog(arange(N),abs(mu))
for n in arange(N):
    mu[n] = sum(f_sample[1,:n+1])
    mu[n] /= (n+1)
pylab.loglog(arange(N),abs(mu))
pylab.show()

# 1D periodic test
'''
i1 = 10
i2 = 50

tp = linspace(0,tp[-1],nt)

N = 5000
f_sample = zeros((2,N))
for n in arange(N):
    tic = time.clock()
    f = simulate.randProcess1DPeriodic(tp, nt, nt)
    toc = time.clock()
    print n,':',toc-tic
    f_sample[0,n] = f[i1]
    f_sample[1,n] = f[i2]

mu1 = sum(f_sample[0,:])/N
mu2 = sum(f_sample[1,:])/N

cov = dot(f_sample[0,:]-mu1,f_sample[1,:]-mu2)/N

print 'estimate:',cov

cov_ex = simulate.covfun(tp[i1]-tp[i2])
print 'exact:', cov_ex

mu = zeros(N)
for n in arange(N):
    mu[n] = sum(f_sample[0,:n+1])
    mu[n] /= (n+1)
pylab.loglog(arange(N),abs(mu))
for n in arange(N):
    mu[n] = sum(f_sample[1,:n+1])
    mu[n] /= (n+1)
pylab.loglog(arange(N),abs(mu))
pylab.show()
'''

'''
i1 = 10
j1 = 10
i2 = 20
j2 = 20

N = 300
cov = 0.0
for n in arange(N):
    tic = time.clock()
    fp = simulate.randProcessPeriodic(sp, tp, ns, nt, Ks, Kt)
    toc = time.clock()
    print n,':',toc-tic
    cov += fp[i1,j1]*fp[i2,j2]

cov /= N

print 'estimate:',cov

th1 = sp[i1]/sp[-1]*2*pi
th2 = sp[i2]/sp[-1]*2*pi
d = sqrt((cos(th1)-cos(th2))**2 + (sin(th1)-sin(th2))**2)
cov_ex = simulate.covfun(d)*simulate.covfun(tp[j1]-tp[j2])
print 'exact:', cov_ex
'''

