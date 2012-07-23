import pylab, os, math, string
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import read_blade
from scipy.interpolate import splprep, splev, sproot, splrep

def sorter(dist, M):    # generate the a list of indices whose distance to isim is nearest to M
    indRecord = []
    sortDist = sorted(dist)
    for i in arange(M+1):
        ind = dist.index(sortDist[i])
        dist[ind] = 999999999
        indRecord += [ind]
    indRecord.pop(0)      # remove the index of itself

    return indRecord

def read_utcfd(rundir):
    # inputs
    # rundir : where runs are stored

    # check if completed
    iran = zeros((len(os.listdir(rundir))))
    i = 0
    for rd in os.listdir(rundir):
        if os.path.exists(rundir+str(rd)+'/utcfd_out.monitoring.dat'):
            lines = file(rundir+str(rd)+'/utcfd_out.monitoring.dat').readlines()
            if lines[-1][1:4] == 'Job':
                iran[i] = 1
            lines = file(rundir+str(rd)+'/Z.dat').readlines()
            nkl = len(lines)
        i += 1

    nruns = int(sum(iran))
    eta  = zeros((nruns))
    pr   = zeros((nruns))
    mdot = zeros((nruns))
    dbeta = zeros((nruns))
    dth  = zeros((nruns))
    Z0   = zeros((nruns,nkl))

    ir = 0
    i = 0
    for rd in os.listdir(rundir):
        if iran[ir]:
            # read in dth (if it exists)
            if os.path.exists(rundir+str(rd)+'/dth.dat'):
                lines = file(rundir+str(rd)+'/dth.dat').readlines()
                dth[i] = array([line.strip().split() for line in lines], float)
            # read in Z values
            lines = file(rundir+str(rd)+'/Z.dat').readlines()
            tmp = array([line.strip().split() for line in lines], float)
            Z0[i,:] = tmp.squeeze()
            # read in inlet quantities
            lines = file(rundir+str(rd)+'/output/utcfd_out.inout_1000.dat').readlines()
            in_vals = array(lines[-1].strip().split(), float)
            # read in outlet quantites
            lines = file(rundir+str(rd)+'/output/utcfd_out.inout_2000.dat').readlines()
            out_vals = array(lines[-1].strip().split(), float)
            # calculate quantities of interest
            mdot[i] = out_vals[1]
            pr[i] = out_vals[2]/in_vals[2]
            tr = out_vals[3]/in_vals[3]
            eta[i] = (pr[i]**(0.4/1.4) - 1.0) / (tr - 1.0)
            dbeta[i] = in_vals[9]-out_vals[9]
            i += 1
        ir += 1

    return eta, pr, mdot, dbeta, Z0, dth

# read in the results
rundir = '/mnt/pwfiles/ericdow/random_blade/runs/'
eta, pr, mdot, dbeta, Z0, dth = read_utcfd(rundir)

# read in the eigenvalues of the K-L
lines = file('/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/kl_modes/S.dat').readlines()
S_kl = array([line.strip().split() for line in lines], float)

# number of simulations to include in stacked Jacobian
Nsim = Z0.shape[0]
NZ = Z0.shape[1]

print 'Total Simulations Available: ', Z0.shape[0]

# K-L mode amplitudes
Z0 = Z0[:Nsim,:NZ]
S_kl = S_kl[:NZ]

# compute the slope of the speedline at design
# mdot1 = 44.408237
# pr1   = 2.06232347
# mdot2 = 44.607510
# pr2   = 2.05032421
mdot1 = 44.456455
pr1   = 2.05916824
mdot2 = 44.557335
pr2   = 2.05348317
alpha = (pr2-pr1)/2.0/(mdot2-mdot1)
print alpha

# matrix of quantities of interest
C1 = eta
C2 = pr + alpha*mdot
# C  = vstack([C1, C2]).T
C = C2

# search for nearest samples
M = 100
for isim in range(Nsim):
    print 'Sorting :', isim
    dist_i = zeros((Nsim))
    for jsim in range(Nsim):         # isim itself is also involved up to now
         # distance**2 between sample isim & jsim, measured in mode space
        dist_i[jsim]  = sum( ( (Z0[isim,:]-Z0[jsim,:]) * S_kl )**2 )    
    # M_index = sorter(dist_i, M)         # the indices of the closest M to isim sample
    M_index = argsort(dist_i)[1:M+1]
    if isim == 0:
        Neighbors = M_index
    else:
        Neighbors = vstack([Neighbors, M_index])    # sample-by-M array

# compute the Jacobian at each sample
for isim in range(Nsim):
    # least squares for each mode
    if C.ndim > 1:
        Y = vstack((C[isim,:],C[Neighbors[isim,:],:]))
    else:
        Y = hstack((C[isim],C[Neighbors[isim,:]]))
    A = ones((M+1,NZ+1))
    A[:,1:] = vstack((Z0[isim,:],Z0[Neighbors[isim,:],:]))
    J = linalg.lstsq(A,Y)[0][1:]
    # form the stacked Jacobian
    if isim == 0:
        stackJ = J.T
    else:
        stackJ = vstack([stackJ, J.T])

# plot average Jacobian for each K-L mode
J_ave = mean(abs(stackJ),axis=0)

pylab.figure()
pylab.semilogy(arange(NZ)+1,J_ave,'*')
pylab.grid(True)
pylab.xlabel('K-L mode index (i)')
pylab.ylabel('Average |dJ/dZ_i|')
pylab.title('Jacobian of speed-line shift with '+str(M)+' nearest neighbors')

# compute SVD of stacked Jacobian
V,S,UT = linalg.svd(stackJ)

# read in the surface mesh
rpath = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/blade_surf.dat'
cdim, sdim, x, y, z = read_blade.read_coords(rpath)

# read in the K-L modes
kl_path = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/kl_modes/'
V_kl = zeros((NZ,cdim,sdim))
for i in range(NZ):
    cdim,sdim,V_kl[i,:,:] = read_blade.read_mode(kl_path+'V'+str(i+1)+'.dat')

V_kl_tmp = zeros((NZ,cdim*sdim))
for i in range(NZ):
    V_kl_tmp[i,:] = reshape(V_kl[i,:,:],(cdim*sdim))

# transform the modes
V_trans_tmp = dot(dot(UT,eye(NZ)*S_kl),V_kl_tmp)

V_trans = zeros((NZ,cdim,sdim))
for i in range(NZ):
    V_trans[i,:,:] = reshape(V_trans_tmp[i,:],(cdim,sdim))

'''
# write out the transformed modes
for imode in range(V_trans.shape[0]):
    fname = 'V'+str(imode+1)+'.dat'
    f = open('transformed_modes/'+fname,'w')
    f.write('CDIM:       %d\n' % cdim)
    f.write('SDIM:       %d\n' % sdim)
    for i in arange(cdim):
        for j in arange(sdim):
            f.write('%e\n' % V_trans[imode,i,j])
    f.close()
'''

# write out the transformed amplitudes
Cov_trans = dot(UT,dot(S_kl**2*eye(NZ),UT.T))
f = open('transformed_modes/S.dat','w')
for i in arange(NZ):
    f.write('%e\n' % sqrt(Cov_trans[i,i]))
f.close()

s, t = read_blade.xyz2st(x,y,z)
s /= s[-1]

# plot dth
pylab.figure()
pylab.hist(dth,bins=50)
pylab.xlabel('dth')

Zs = tile(t,(cdim,1))
# chordwise location of slices
Cs = zeros((sdim,cdim))
for isec in arange(sdim):
    tck,sn = splprep([x[:,isec],y[:,isec],z[:,isec]],s=0,per=1)
    Cs[isec,:] = sn

Cs = 2.0*(Cs - 0.5)
CsT = (1./pi*arccos(1.0-2.0*(-abs(Cs)+1.0)) - 1.0)*0.5*(1.0-sign(Cs))\
     + 1./pi*arccos(1.0-2.0*abs(Cs))*0.5*(1.0+sign(Cs))
CsT = 0.5*CsT + 0.5
Cs  = 0.5*Cs  + 0.5

# plot the transformed modes
# make sure they match up with the K-L modes
ls = len(s)+mod(len(s),2)
V_trans_tmp = zeros(V_trans.shape)
V_trans_tmp[:,:ls/2,:] = V_trans[:,ls/2-1:,:]
V_trans_tmp[:,ls/2-1:,:] = V_trans[:,:ls/2,:]
for i in range(4):
    pylab.figure()
    pylab.contourf(s,t,V_trans_tmp[i,:,:].T,50)
    # pylab.contourf(CsT,Zs.T,V_trans_tmp[i,:,:].T,50)
    pylab.title('Output-based mode '+str(i+1))
    pylab.xlabel('Chordwise location')
    pylab.ylabel('Spanwise Coordinate')
    pylab.text(0.03,t[0]-0.07,'LE')
    pylab.text(0.5,t[0]-0.07,'TE')
    pylab.text(0.97,t[0]-0.07,'LE')
    pylab.text(0.25,t[0]-0.07,'PS')
    pylab.text(0.75,t[0]-0.07,'SS')
    pylab.ylim([t[0]-0.1,t[-1]])
    pylab.colorbar()

# plot the decay of the eigenvalues of the stacked Jacobian
pylab.figure()
pylab.semilogy(arange(NZ)+1,S**2,'*')
pylab.title('Eigenvalues of output-based modes')
pylab.grid('True')
pylab.xlim((0,NZ+1))
pylab.xlabel('Mode Index')
pylab.ylabel('Eigenvalue')

pylab.figure()
pylab.semilogy(arange(NZ)+1,S_kl,'*')
pylab.title('Eigenvalues of K-L Modes')
pylab.grid('True')
pylab.xlim((0,NZ+1))
pylab.xlabel('Mode Index')
pylab.ylabel('Eigenvalue')

# plot the eifficiency against the first mode amplitude
W = zeros((Nsim,NZ))
for i in range(Nsim):
    W[i,:] = dot(dot(Z0[i,:],S_kl*eye(NZ)),UT.T)

# use regression to evaluate nonlinearity
A = ones((Nsim,6))
Y = copy(eta)
A[:,1] = W[:,0]
A[:,2] = W[:,1]
A[:,3] = W[:,0]*W[:,1]
A[:,4] = W[:,0]**2
A[:,5] = W[:,1]**2
c = linalg.lstsq(A,Y)[0]

# plot the regression surface
W1 = linspace(W[:,0].min(),W[:,0].max(),100)
W2 = linspace(W[:,1].min(),W[:,1].max(),100)
eta_reg = zeros((100,100))
for i in range(100):
    for j in range(100):
        eta_reg[i,j] = c[0] + c[1]*W1[i] + c[2]*W2[j] + c[3]*W1[i]*W2[j] + c[4]*W1[i]**2 + c[5]*W2[j]**2
pylab.figure()
pylab.contourf(W1,W2,eta_reg.T,50)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Efficiency Regression Surface')

pylab.figure()
img = pylab.contourf(W1,W2,eta_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W[:,0],W[:,1],100,eta,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.title('Efficiency Regression Surface')
    
pylab.figure()
pylab.plot(W[:,0],eta,'*')
pylab.grid(True)
pylab.xlabel('Dominant Mode Amplitude')
pylab.ylabel('Efficiency')
pylab.title('Efficiency vs Dominant Mode Amplitude')
pylab.savefig('figures/eta_vs_dom.png')

pylab.figure()
pylab.plot(W[:,0],pr,'*')
pylab.grid(True)
pylab.xlabel('Dominant Mode Amplitude')
pylab.ylabel('Pressure Ratio')
pylab.title('Pressure Ratio vs Dominant Mode Amplitude')
pylab.savefig('figures/pr_vs_dom.png')

pylab.figure()
pylab.plot(W[:,0],mdot,'*')
pylab.grid(True)
pylab.xlabel('Dominant Mode Amplitude')
pylab.ylabel('Mass Flow Rate [lb/s]')
pylab.title('Mass Flow Rate vs Dominant Mode Amplitude')
pylab.savefig('figures/mdot_vs_dom.png')

pylab.figure()
pylab.plot(W[:,0],dbeta,'*')
pylab.grid(True)
pylab.xlabel('Dominant Mode Amplitude')
pylab.ylabel('Turning Angle [degrees]')
pylab.title('Turning Angle vs Dominant Mode Amplitude')
pylab.savefig('figures/dbeta_vs_dom.png')

pylab.figure()
pylab.plot(dbeta,eta,'*')
pylab.grid(True)
pylab.xlabel('Turning Angle [degrees]')
pylab.ylabel('Efficiency')
pylab.title('Turning Angle vs Efficiency')
pylab.savefig('figures/dbeta_vs_eta.png')

pylab.figure()
pylab.plot(pr,eta,'*')
pylab.grid(True)
pylab.xlabel('Pressure Ratio')
pylab.ylabel('Efficiency')
pylab.title('Pressure Ratio vs Efficiency')
pylab.savefig('figures/pr_vs_eta.png')

pylab.figure()
pylab.plot(mdot,eta,'*')
pylab.grid(True)
pylab.xlabel('Mass Flow Rate [lb/s]')
pylab.ylabel('Efficiency')
pylab.title('Mass Flow Rate vs Efficiency')
pylab.savefig('figures/mdot_vs_eta.png')

pylab.show()

