import pylab
from numpy import *
from numpy.linalg import *
from scipy import stats

N = 100
th = linspace(0,2*pi,N)

x = cos(th)
y = sin(th)
d = sqrt((x - x[0])**2 + (y - y[0])**2)

sig = 1.0
lam = 0.25
f = sig**2*exp(-d**2/(2*lam**2))

F = fft.fft(f).real/N

Z = random.randn(N)
p = zeros((1,N))
for i in arange(N/2+1):
    p += Z[i]*F[i]*cos(i*th)
for i in arange(1,N/2):
    p += Z[N/2+i]*F[i]*sin(i*th)

pylab.plot(th,p.T)
pylab.plot(th+2*pi,p.T)
pylab.plot(th[-1],p.T[-1],'ro')
pylab.show()


