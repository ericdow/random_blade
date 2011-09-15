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
    lam, sig = 1.0, sqrt(1.0e-2)
    return (sig**2) * exp(-x**2/(2*lam**2))

def cov(x, covF): 
    covM = array([[covF(a-b) for a in x] for b in x])
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
    N = covM.shape[0]
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

def randProcess1D(y, ny, KY):   
    muY = mu(y, mufun)
    CY  = cov(y, covfun)
    UY, SY = eigenPairs(CY, y[0], y[-1])

    '''
    # check that these are eigenfunctions
    check = 5
    U_int = zeros((ny))
    dy = abs((y[-1]-y[0])/ny)
    for i in arange(len(y)):
        Ker = array([covfun(y[i]-a) for a in y])
        U_int[i] = dy*dot(UY[:,check],Ker)

    pylab.plot(y,U_int)
    pylab.plot(y,SY[check]*UY[:,check],'o')
    print (SY[check]*UY[:,check]) / U_int
    pylab.show()
    '''
    
    random.seed() 
    Z = random.randn(KY)    # normal distribution vector
    Z *= sqrt(SY[:KY])
    f = dot(UY[:,:KY],Z)
    f += muY
    return f

def randProcess1DPeriodic(x, nx, KX):

    '''
    x -= x[0]
    th = linspace(0,2*pi,nx)
    xx = x[-1]*cos(th)/2/pi # scale by radius 
    yy = x[-1]*sin(th)/2/pi
    d = sqrt((xx - xx[0])**2 + (yy - yy[0])**2)
    F = fft.fft(covfun(d)).real/nx
    for i in arange(len(F)):
        if F[i] < 0.0:
            F[i:] = 0.0
            break

    # F = F/nx*x[-1]
    
    random.seed() 
    Z = random.randn(nx)    # normal distribution vector
    f = zeros(nx)
    for i in arange(nx/2+1):
        f += sqrt(F[i])*Z[i]*cos(i*th*2*pi/x[-1])*sqrt(2/x[-1])
    for i in arange(1,nx/2):
        f += sqrt(F[i])*Z[nx/2+i]*sin(i*th*2*pi/x[-1])*sqrt(2/x[-1])

    
    # check that these are eigenfunctions
    check = 2
    U_int = zeros((nx))
    dx = abs(x[-1]/nx)
    ef1 = cos(check*x*2*pi/x[-1])*sqrt(2/x[-1])
    ef2 = sin(check*x*2*pi/x[-1])*sqrt(2/x[-1])
    for i in arange(len(x)):
        Ker = array([covfun(sqrt((cos(2*pi*x[i]/x[-1]) - cos(2*pi*a/x[-1]))**2\
                + (sin(2*pi*x[i]/x[-1]) - sin(2*pi*a/x[-1]))**2)) for a in x])
        U_int[i] = dx*dot(ef2,Ker)


    pylab.plot(x,U_int)
    pylab.plot(x,F[check]*ef2,'r--')
    print (F[check]*ef2) / U_int
    pylab.show()
    '''
    muX = mu(x, mufun)
    CX = covPeriodic(x,covfun)
    UX, SX = eigenPairs(CX, x[0], x[-1])

    '''
    # check that these are eigenfunctions
    check = 0
    U_int = zeros((nx))
    dx = abs((x[-1]-x[0])/nx)
    for i in arange(len(x)):
        Ker = array([covfun(sqrt((cos(2*pi*x[i]/x[-1]) - cos(2*pi*a/x[-1]))**2\
                + (sin(2*pi*x[i]/x[-1]) - sin(2*pi*a/x[-1]))**2)) for a in x])
        U_int[i] = dx*dot(UX[:,check],Ker)

    pylab.plot(x,U_int)
    pylab.plot(x,SX[check]*UX[:,check],'o')
    print (SX[check]*UX[:,check]) / U_int
    pylab.show()
    '''

    random.seed() 
    Z = random.randn(KX)    # normal distribution vector
    Z *= sqrt(SX[:KX])
    f = dot(UX[:,:KX],Z)
    f += muX
    
    return f

'''
def randProcessPeriodic(x, y, nx, ny, KX, KY):
    # process is periodic in the x direction
    muX = mu(x, mufun)
    CX = covPeriodic(x,covfun)
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
'''

def randProcessPeriodic(x, y, nx, ny, KX, KY):
    # process is periodic in the x direction
    xl = linspace(0.,x[-1],nx)
    yl = linspace(0.,y[-1],ny)
    muX = mu(xl, mufun)
    CX = covPeriodic(xl,covfun)
    UX, SX = eigenPairs(CX, xl[0], xl[-1])
    # pylab.figure()
    # pylab.plot(xl,UX[:,0])
    for i in arange(KX):
        u = interp1d(xl, UX[:,i])
        UX[:,i] = u(x)
    # pylab.plot(x,UX[:,0])
    
    muY = mu(yl, mufun)
    CY  = cov(yl, covfun)
    UY, SY = eigenPairs(CY, yl[0], yl[-1])
    # pylab.figure()
    # pylab.plot(yl,UY[:,0])
    for i in arange(KY):
        u = interp1d(yl, UY[:,i])
        UY[:,i] = u(y)
    # pylab.plot(y,UY[:,0])

    # pylab.show()

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
            
    '''
    nx = len(x)
    ny = len(y)
    nxc = len(xc)
    nyc = len(yc)
    xg, yg = meshgrid(x,y)
    xcg, ycg = meshgrid(xc,yc)
    xg = xg.reshape(nx*ny,1)
    yg = yg.reshape(nx*ny,1)
    xcg = xcg.reshape(nxc*nyc,1)
    ycg = ycg.reshape(nxc*nyc,1)
    f = f.reshape(nx*ny,1)
    # fint = interp2d(xg,yg,f)
    fint = interp2d(x,y,f)
    ff = zeros((len(xcg),1))
    for i in arange(nxc):
        print i
        ff[i] = fint(xcg[i],ycg[i])
    ff = ff.reshape(nxc,nyc)
    return ff
    '''



 
