from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import pylab
from numpy import *
from numpy.linalg import *
from scipy import stats
import math
from scipy.interpolate import interp2d

def mufun(x,y):    # expectation function
    f = 0
    return f

def covfun(dx,dy):    # covariance function
    lam, sig = 0.2, 1e-4
    return (sig**2) * exp( -(dx**2 + dy**2)/ (lam**2) )

def covMatrix(covF, x, y):    # generate covariance matrix on points x
    covM = array([[covF(x[i]-x[j],y[i]-y[j]) for i in arange(len(x))]\
                                             for j in arange(len(x))])
    return covM

def covGen(x, y, x0, y0, covF): 
#   x,y:    all points being studied
#   x0,y0:  measured points
#   covF:   cov function handle
#   Cc:     returned conditioned covariance matrix
    C = covMatrix(covF, x, y).squeeze()
    b = array([[covF(x0[i]-x[j],y0[i]-y[j]) for i in arange(len(x0))]\
                                            for j in arange(len(x))])
    b = b.squeeze()
    D = covMatrix(covF, x0, y0).squeeze()
    Cc= C - dot(b, linalg.solve(D,b.T)) # dot(linalg.inv(D), b.T))
    return Cc

def muGen(x, y, x0, y0, d, muF, covF):
#   x,y:    all points being studied
#   x0,y0:  measured points
#   d:      values measured at (x0, y0)
#   muF:    mu function handle
#   covF:   cov function handle
#   muC:    returned conditioned mu(x)
    b = array([[covF(x0[i]-x[j],y0[i]-y[j]) for i in arange(len(x0))]\
                                            for j in arange(len(x))])
    b = b.squeeze()
    D = covMatrix(covF, x0, y0).squeeze()
    muC = dot(b, linalg.solve(D,(d-muF(x0,y0)).T))
    muC += muF(x,y)
    return muC 

def eigenPairs(covM, a, b):    # generate eigen-pairs, x range is [a,b]
    U, S = linalg.svd(covM)[:2]
    N = covM.shape[0]
    U = U * sqrt(  N / (b - a)  );
    S = S / N;
    return U, S

def randProcess(mu, U, S, K, a, b):    # generate unconditioned random process f(x,w) using K-L expansion
#   mu: vector
#   K:  high mode truncation, using modes 0,1,2,...K-1
#   f:  returned random process
    random.seed()    # initialize random number generator
    Z = random.randn(K)    # normal distribution vector
    Z *= sqrt(S[:K])
    f = dot(U[:,:K],Z)
    f += mu

    return f

def grid(nx, ny, xmax, ymax):
    x, y = meshgrid(linspace(0,xmax,nx),linspace(0,ymax,ny))
    x = x.reshape(nx*ny,1)
    y = y.reshape(nx*ny,1)

    return x, y

def interp(x, y, f, xc, yc):
    fint = interp2d(x,y,f)
    ff = zeros((len(xc),1))
    for i in arange(len(xc)):
        ff[i] = fint(xc[i],yc[i])
    
    return ff
        
xmax = 3.
ymax = xmax
ny = 60
nx = ny
x, y = grid(nx, ny, xmax, ymax)

a = 0
b = xmax
K = nx*ny/5          # mode cutoff

y0 = hstack((linspace(0,ymax,ny),linspace(0,ymax,ny))).T
x0 = vstack((zeros((ny,1)),xmax*ones((ny,1))))
d  = zeros(x0.size)

muC  =  muGen(x, y, x0, y0, d, mufun, covfun)
Cc   = covGen(x, y, x0, y0, covfun)
U, S = eigenPairs(Cc, a, b)
f = randProcess(muC, U, S, K, a, b)

x = x.reshape(ny,nx)
y = y.reshape(ny,nx)
f = f.reshape(ny,nx)

xg = x[0,:]
yg = y[:,0]
xint = array([0., 0.5, 1., 1.5])
yint = array([0., 0.5, 1., 1.5])
fint = interp(xg, yg, f, xint, yint)

fig = pylab.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(x, y, f, rstride=1, cstride=1, cmap = cm.jet)
# fig.colorbar(surf)
pylab.show()
