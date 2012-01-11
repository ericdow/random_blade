from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
from numpy import *
from numpy.linalg import *
from scipy import stats
import math
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d

def mufun(x):    # expectation function
    f = 0
    return f

def covfun(x):    # covariance function
    lam, sig = 2.0, sqrt(20e-3)
    return (sig**2) * exp(-x**2/(2*lam**2))

def covfun_span(x,s_params):    # covariance function
    lam = s_params[0]
    sig = s_params[1]
    return (sig**2) * exp(-x**2/(2*lam**2))

def covfun_chord(x1,x2,c_params):    # covariance function
    # stationary parameters
    lam = c_params[0]
    # weighting parameters
    sig = c_params[1]
    w = c_params[2]
    s = c_params[3]
    xc = c_params[4]

    if x1 < xc - w/2.0:
        a1 = sig*exp(-(x1-(xc-w/2.0))**2/s**2)
    elif x1 > xc + w/2.0:
        a1 = sig*exp(-(x1-(xc+w/2.0))**2/s**2)
    else:
        a1 = sig

    if x2 < xc - w/2.0:
        a2 = sig*exp(-(x2-(xc-w/2.0))**2/s**2)
    elif x2 > xc + w/2.0:
        a2 = sig*exp(-(x2-(xc+w/2.0))**2/s**2)
    else:
        a2 = sig
        
    return a1 * a2 * exp(-(x1-x2)**2/(2*lam**2))

def covfun2(x):    # covariance function
    lam, sig = 2.0, sqrt(20e-3)
    return (sig**2) * exp(-x**2/(2*lam**2))

# stationary covariance
def cov(x, covF, *params): 
    covM = array([[covF(a-b,*params) for a in x] for b in x])
    return covM

# nonstationary covariance
def covNonstat(x, covF, *params): 
    covM = array([[covF(a,b,*params) for a in x] for b in x])
    return covM

def covPeriodic(x, covF):
    xt = 2/pi*x
    covM = array([[covfun(sqrt((cos(2*pi*a/xt[-1])-cos(2*pi*b/xt[-1]))**2\
                           + (sin(2*pi*a/xt[-1])-sin(2*pi*b/xt[-1]))**2))\
                   for a in xt] for b in xt])
    return covM

def covCond(x, x0, covF): 
    C = cov(x, covF)    # unconditioned covariance matrix
    b = array([[covF(a-b) for a in x0] for b in x])
    D = cov(x0, covF)
    Cc= C - dot(b, linalg.solve(D,b.T)) 
    return Cc

def mu(x, muF):
    return array([muF(a) for a in x])

def muCond(x, x0, d, muF, covF):
    b = array([[covF(a-b) for a in x0] for b in x])
    D = cov(x0, covF)
    muC = dot(b, linalg.solve(D,(d-mu(x0,muF)).T))
    muC += mu(x, muF)
    return muC 

def eigenPairs(covM, a, b):    # generate eigen-pairs, x range is [a,b]
    U, S = linalg.svd(covM)[:2]
    # make sure unique sign
    N = covM.shape[0]
    for i in arange(N):
        U[:,i] *= sign(U[1,i])
    U = U * sqrt(  N / (b - a)  );
    S = S / N * (b-a);
    return U, S

def randProcessCond(x, y, x0, d, nx, ny, KX, KY):
    muX = muCond(x, x0, d, mufun, covfun)
    CX  = covCond(x, x0, covfun)
    UX, SX = eigenPairs(CX, x[0], x[-1])
    
    muY = mu(y, mufun)
    CY  = cov(y, covfun)
    UY, SY = eigenPairs(CY, y[0], y[-1])

    random.seed() 
    Z = random.randn(KX,KY)    # normal distribution vector
    f = zeros((nx,ny))
    for i in arange(KX):
        for j in arange(KY):
            f += sqrt(SX[i])*sqrt(SY[j])*Z[i,j]*outer(UX[:,i],UY[:,j])

    f += outer(muX,ones((1,ny)))
    f += outer(muY,ones((1,nx))).T
    return f

def randProcess(x, y, nx, ny, KX, KY):
    muX = mu(x, mufun)
    CX  = cov(x, covfun)
    UX, SX = eigenPairs(CX, x[0], x[-1])
    
    muY = mu(y, mufun)
    CY  = cov(y, covfun)
    UY, SY = eigenPairs(CY, y[0], y[-1])

    # REMOVE SEED NUMBER TO MAKE RANDOM
    random.seed(5) 
    Z = random.randn(KX,KY)    # normal distribution vector
    f = zeros((nx,ny))
    for i in arange(KX):
        for j in arange(KY):
            f += sqrt(SX[i])*sqrt(SY[j])*Z[i,j]*outer(UX[:,i],UY[:,j])
    
    f += outer(muX,ones((1,ny)))
    f += outer(muY,ones((1,nx))).T
    return f

def randProcessPeriodic(x, y, nx, ny, KX, KY):
    # process is periodic in the x direction
    xl = linspace(0.,x[-1],nx)
    yl = linspace(0.,y[-1],ny)
    muX = mu(xl, mufun)
    CX = covPeriodic(xl,covfun)
    UX, SX = eigenPairs(CX, xl[0], xl[-1])
    # interpolate from uniform to original grid
    for i in arange(KX):
        u = interp1d(xl, UX[:,i])
        UX[:,i] = u(x)
    
    muY = mu(yl, mufun)
    CY  = cov(yl, covfun2)
    UY, SY = eigenPairs(CY, yl[0], yl[-1])
    # interpolate from uniform to original grid
    for i in arange(KY):
        u = interp1d(yl, UY[:,i])
        UY[:,i] = u(y)

    # REMOVE SEED NUMBER TO MAKE RANDOM
    random.seed(7) 
    Z = random.randn(KX,KY)    # normal distribution vector
    f = zeros((nx,ny))
    for i in arange(KX):
        for j in arange(KY):
            f += sqrt(SX[i])*sqrt(SY[j])*Z[i,j]*outer(UX[:,i],UY[:,j])

    f += outer(muX,ones((1,ny)))
    f += outer(muY,ones((1,nx))).T
    return f
    
def randProcessLE(x, y, i_le, nx, ny, KX, KY, c_params, s_params):
    # x direction is chord
    # y direction is span
    # i_le gives the index of the LE coordinate
    muX = mu(x, mufun)
    CX  = covNonstat(x, covfun_chord, c_params)
    UX, SX = eigenPairs(CX, x[0], x[-1])

    # plot the eigenfunctions
    '''
    pylab.figure()
    for i in arange(5):
        pylab.plot(x,UX[:,i])
    pylab.show()
    '''
    
    muY = mu(y, mufun)
    CY  = cov(y, covfun_span, s_params)
    UY, SY = eigenPairs(CY, y[0], y[-1])

    random.seed() 
    Z = random.randn(KX,KY)    # normal distribution vector
    f = zeros((nx,ny))
    for i in arange(KX):
        for j in arange(KY):
            f += sqrt(SX[i])*sqrt(SY[j])*Z[i,j]*outer(UX[:,i],UY[:,j])
    
    f += outer(muX,ones((1,ny)))
    f += outer(muY,ones((1,nx))).T
    return f

def grid(nx, ny, xmax, ymax):
    x, y = meshgrid(linspace(0,ymax,ny),linspace(0,xmax,nx))
    return x, y

def interp(x, y, f, xc, yc):
    # interpolate f from x/y to xc/yc
    ff = zeros((len(xc),len(yc)))
    for i in arange(len(xc)):
        for j in arange(len(yc)):
            # find nearest NxN points
            N = 6
            xm = nonzero(x<xc[i])
            xm = xm[0]
            if len(xm) == 0:
                xis = 0
                xie = N
            elif xm[-1] < (N/2-1):
                xis = 0
                xie = N
            elif xm[-1] > len(x)-(N/2+1):
                xis = len(x)-(N+1)
                xie = len(x)-1
            else:
                xis = xm[-1]-(N/2-1)
                xie = xm[-1]+(N/2+1)
                
            ym = nonzero(y<yc[j])
            ym = ym[0]
            if len(ym) == 0:
                yis = 0
                yie = N
            elif ym[-1] < (N/2-1):
                yis = 0
                yie = N
            elif ym[-1] > len(y)-(N/2+1):
                yis = len(y)-(N+1)
                yie = len(y)-1
            else:
                yis = ym[-1]-(N/2-1)
                yie = ym[-1]+(N/2+1)

            # print xis,xie,yis,yie
            # print f.shape
            fi = f[xis:xie,yis:yie]
            fint = interp2d(x[xis:xie],y[yis:yie],fi,kind='linear')
            ff[i,j] = fint(xc[i],yc[j])
            ff[i,j] = f[(xis+xie)/2,(yis+yie)/2]

    return ff
 
