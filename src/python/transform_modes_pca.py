import pylab, os, math, string
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.stats.stats import pearsonr
import read_blade, write_tecplot

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
            npca = len(lines)
        i += 1

    nruns = int(sum(iran))
    eta  = zeros((nruns))
    pr   = zeros((nruns))
    mdot = zeros((nruns))
    dbeta = zeros((nruns))
    Z0   = zeros((nruns,npca))

    ir = 0
    i = 0
    for rd in os.listdir(rundir):
        if iran[ir]:
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

    return eta, pr, mdot, dbeta, Z0

def read_utcfd_meas(rundir):
    # inputs
    # rundir : where runs are stored

    # read in inlet quantities
    lines = file(rundir+'/output/utcfd_out.inout_1000.dat').readlines()
    in_vals = array(lines[-1].strip().split(), float)
    # read in outlet quantites
    lines = file(rundir+'/output/utcfd_out.inout_2000.dat').readlines()
    out_vals = array(lines[-1].strip().split(), float)
    # calculate quantities of interest
    mdot = out_vals[1]
    pr = out_vals[2]/in_vals[2]
    tr = out_vals[3]/in_vals[3]
    eta = (pr**(0.4/1.4) - 1.0) / (tr - 1.0)
    dbeta = in_vals[9]-out_vals[9]

    return eta, pr, mdot, dbeta

# compute the slope of the speedline at design
mdot1 = 44.456455
pr1   = 2.05916824
mdot2 = 44.557335
pr2   = 2.05348317
alpha = (pr2-pr1)/2.0/(mdot2-mdot1)

# read in the nominal performance    
nomdir = '/mnt/pwfiles/ericdow/run/coarse/stationary_hub/'
eta_nom, pr_nom, mdot_nom, dbeta_nom = read_utcfd_meas(nomdir)
sls_nom = pr_nom - alpha*mdot_nom

# read in the results
rundir = '/mnt/pwfiles/ericdow/random_blade/runs/'
eta, pr, mdot, dbeta, Z0 = read_utcfd(rundir)

'''
# compute operating point shift in tangent direction vs orthogonal direction
v1 = array((1,alpha))
v1 = v1/linalg.norm(v1) # unit vector tangent to speedline
v2 = zeros(v1.shape)
v2[0] = -v1[1]
v2[1] = v1[0] # unit vector orthogonal to speedline
d = vstack((mdot-mdot_nom,pr-pr_nom)).T
d1 = dot(d,v1) # shift in the v1 direction
d2 = dot(d,v2) # shift in the v2 direction
pylab.plot(d1,d2,'x')
pylab.axis('Equal')
pylab.grid(True)
pylab.xlabel('Shift in tangent direction')
pylab.ylabel('Shift in orthogonal direction')

pylab.figure()
pylab.plot(eta,pr,'x')
pylab.plot(eta_nom,pr_nom,'ko')
pylab.xlabel('Efficiency')
pylab.ylabel('Pressure ratio')
pylab.grid(True)
'''

# read in the eigenvalues of the PCA
lines = file('/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/pca_modes/S.dat').readlines()
S_pca = array([line.strip().split() for line in lines], float)

# number of simulations to include in stacked Jacobian
Nsim = Z0.shape[0]
NZ = Z0.shape[1]

print 'Total Simulations Available: ', Z0.shape[0]

# PCA mode amplitudes
Z0 = Z0[:Nsim,:NZ]
S_pca = S_pca[:NZ]

# matrix of quantities of interest
C1 = eta
C2 = pr - alpha*mdot
C3 = mdot
C1_nom = eta_nom
C2_nom = sls_nom
C3_nom = mdot_nom
C  = vstack([C1, C2]).T
C = C3

'''
dm = (mdot-mdot_nom)/mdot_nom
dm = (pr - pr_nom)/pr_nom
dS = (C2-sls_nom)/sls_nom
pylab.plot(dm,dS,'x')
pylab.grid(True)
pylab.axis('Equal')
pylab.xlabel('dm/m0')
pylab.ylabel('dS/S0')

r = pearsonr(dm,dS)
'''

# search for nearest samples
M = 100
for isim in range(Nsim):
    print 'Sorting :', isim
    dist_i = zeros((Nsim))
    for jsim in range(Nsim):         # isim itself is also involved up to now
         # distance**2 between sample isim & jsim, measured in mode space
        dist_i[jsim]  = sum( ( (Z0[isim,:]-Z0[jsim,:]) * S_pca )**2 )    
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

# mean of the Jacobians
meanJ = mean(stackJ,axis=0)

# compute SVD of stacked Jacobian
V,S,UT = linalg.svd(stackJ)

# read in the surface mesh
rpath = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/blade_surf.dat'
cdim, sdim, x, y, z = read_blade.read_coords(rpath)

# read in the PCA modes
pca_path = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/pca_modes_gp/'
V_pca = zeros((NZ,cdim,sdim))
for i in range(NZ):
    cdim,sdim,V_pca[i,:,:] = read_blade.read_mode(pca_path+'V'+str(i+1)+'.dat')

V_pca_tmp = zeros((NZ,cdim*sdim))
for i in range(NZ):
    V_pca_tmp[i,:] = reshape(V_pca[i,:,:],(cdim*sdim))

# transform the modes
V_trans_tmp = dot(dot(UT.T,eye(NZ)*S_pca),V_pca_tmp)

V_trans = zeros((NZ,cdim,sdim))
for i in range(NZ):
    V_trans[i,:,:] = reshape(V_trans_tmp[i,:],(cdim,sdim))

# write out the transformed modes
for imode in range(V_trans.shape[0]):
    fname = 'V'+str(imode+1)+'.dat'
    f = open('ob_modes/'+fname,'w')
    f.write('CDIM:       %d\n' % cdim)
    f.write('SDIM:       %d\n' % sdim)
    for i in arange(cdim):
        for j in arange(sdim):
            f.write('%e\n' % V_trans[imode,i,j])
    f.close()

# write out the transformed amplitudes
Cov_trans = dot(UT,dot(S_pca**2*eye(NZ),UT.T))
f = open('ob_modes/S.dat','w')
for i in arange(NZ):
    f.write('%e\n' % sqrt(Cov_trans[i,i]))
f.close()

s, t = read_blade.xyz2st(x,y,z)
s /= s[-1]

# tecplot file with ob modes
write_tecplot.write_blade_surf(x,y,z,V_trans,'ob_modes_r37.dat') 

# plot the transformed modes
# make sure they match up with the PCA modes
ls = len(s)+mod(len(s),2)
V_trans_tmp = zeros(V_trans.shape)
V_trans_tmp[:,:ls/2,:] = V_trans[:,ls/2-1:,:]
V_trans_tmp[:,ls/2-1:,:] = V_trans[:,:ls/2,:]
for i in range(4):
    pylab.figure()
    pylab.contourf(s,t,V_trans_tmp[i,:,:].T,50)
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
pylab.semilogy(arange(NZ)+1,S**2/S[0]**2,'*')
pylab.title('Normalized eigenvalues of Jacobian covariance')
pylab.grid('True')
pylab.xlim((0,NZ+1))
pylab.ylim((0.9*S[-1]**2/S[0]**2,1.1))
pylab.xlabel('Mode Index')
pylab.ylabel('Eigenvalue')

'''
pylab.figure()
pylab.semilogy(arange(NZ)+1,S_pca**2,'*')
pylab.title('Eigenvalues of PCA Modes')
pylab.grid('True')
pylab.xlim((0,NZ+1))
pylab.xlabel('Mode Index')
pylab.ylabel('Eigenvalue')
'''

# plot the eifficiency against the first mode amplitude
# W = Z0^T*sqrt(Lambda_PCA)*U
W = zeros((Nsim,NZ))
for i in range(Nsim):
    W[i,:] = dot(UT,Z0[i,:])

# use regression to evaluate nonlinearity
A = ones((Nsim,6))
Y = copy(C1)
A[:,1] = W[:,0]
A[:,2] = W[:,1]
A[:,3] = W[:,0]*W[:,1]
A[:,4] = W[:,0]**2
A[:,5] = W[:,1]**2
c_C1 = linalg.lstsq(A,Y)[0]

'''
# plot the regression surface
W1 = linspace(W[:,0].min(),W[:,0].max(),100)
W2 = linspace(W[:,1].min(),W[:,1].max(),100)
eta_reg = zeros((100,100))
for i in range(100):
    for j in range(100):
        eta_reg[i,j] = c_C1[0] + c_C1[1]*W1[i] + c_C1[2]*W2[j] +\
                       c_C1[3]*W1[i]*W2[j] + c_C1[4]*W1[i]**2 + c_C1[5]*W2[j]**2

pylab.figure()
img = pylab.contourf(W1,W2,eta_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W[:,0],W[:,1],100,C1,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Efficiency Regression Surface')
'''

A = ones((Nsim,6))
Y = copy(C2)
A[:,1] = W[:,0]
A[:,2] = W[:,1]
A[:,3] = W[:,0]*W[:,1]
A[:,4] = W[:,0]**2
A[:,5] = W[:,1]**2
c_C2 = linalg.lstsq(A,Y)[0]

'''                       
# plot the regression surface
W1 = linspace(W[:,0].min(),W[:,0].max(),100)
W2 = linspace(W[:,1].min(),W[:,1].max(),100)
sls_reg = zeros((100,100))
for i in range(100):
    for j in range(100):
        sls_reg[i,j] = c_C2[0] + c_C2[1]*W1[i] + c_C2[2]*W2[j] +\
                       c_C2[3]*W1[i]*W2[j] + c_C2[4]*W1[i]**2 + c_C2[5]*W2[j]**2

pylab.figure()
img = pylab.contourf(W1,W2,sls_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W[:,0],W[:,1],100,C2,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Speedline Shift Regression Surface')
'''

A = ones((Nsim,6))
Y = copy(C3)
A[:,1] = W[:,0]
A[:,2] = W[:,1]
A[:,3] = W[:,0]*W[:,1]
A[:,4] = W[:,0]**2
A[:,5] = W[:,1]**2
c_C3 = linalg.lstsq(A,Y)[0]

C1_samp_reg = c_C1[0] + c_C1[1]*W[:,0] + c_C1[2]*W[:,1] +\
               c_C1[3]*W[:,0]*W[:,1] + c_C1[4]*W[:,0]**2 + c_C1[5]*W[:,1]**2
C2_samp_reg = c_C2[0] + c_C2[1]*W[:,0] + c_C2[2]*W[:,1] +\
               c_C2[3]*W[:,0]*W[:,1] + c_C2[4]*W[:,0]**2 + c_C2[5]*W[:,1]**2
C3_samp_reg = c_C3[0] + c_C3[1]*W[:,0] + c_C3[2]*W[:,1] +\
               c_C3[3]*W[:,0]*W[:,1] + c_C3[4]*W[:,0]**2 + c_C3[5]*W[:,1]**2

pylab.figure()
pylab.hist(C1-C1_nom,bins=50,facecolor='y')
pylab.hist(C1-C1_samp_reg,bins=5,facecolor='b')
pylab.xlabel('Efficiency')
pylab.ylabel('Frequency')

pylab.figure()
pylab.hist(C2-C2_nom,bins=50,facecolor='y')
pylab.hist(C2-C2_samp_reg,bins=50,facecolor='b')
pylab.xlabel('Speedline shift')
pylab.ylabel('Frequency')

pylab.figure()
pylab.hist(C3-C3_nom,bins=50,facecolor='y')
pylab.hist(C3-C3_samp_reg,bins=50,facecolor='b')
pylab.xlabel('Mass flowrate')
pylab.ylabel('Frequency')

pylab.figure()
xx = linspace(min(C1),max(C1),2)
pylab.plot(xx,xx,'--')
pylab.plot(C1_samp_reg,C1,'kx',mew=2)
pylab.xlabel('Predicted efficiency')
pylab.ylabel('Actual efficiency')
pylab.grid(True)
pylab.xlim((min(C1),max(C1)))
pylab.ylim((min(C1),max(C1)))

pylab.figure()
xx = linspace(min(C2),max(C2),2)
pylab.plot(xx,xx,'--')
pylab.plot(C2_samp_reg,C2,'kx',mew=2)
pylab.xlabel('Predicted speedline shift')
pylab.ylabel('Actual speedline shift')
pylab.grid(True)
pylab.xlim((min(C2),max(C2)))
pylab.ylim((min(C2),max(C2)))

pylab.figure()
xx = linspace(min(C3),max(C3),2)
pylab.plot(xx,xx,'--')
pylab.plot(C3_samp_reg,C3,'kx',mew=2)
pylab.xlabel('Predicted mass flowrate')
pylab.ylabel('Actual mass flowrate')
pylab.grid(True)
pylab.xlim((min(C3),max(C3)))
pylab.ylim((min(C3),max(C3)))

'''
nb = 34
eta_meas   = zeros((nb))
pr_meas    = zeros((nb))
mdot_meas  = zeros((nb))
dbeta_meas = zeros((nb))
Z_meas     = zeros((nb,NZ))
W_meas     = zeros((nb,NZ))
for i in range(nb):
    # read in the Z values for each blade
    lines = file('/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/pca_modes/Z_blade'+str(i+1)+'.dat').readlines()
    Z_meas[i,:] = array([line.strip().split() for line in lines], float).squeeze()
    # compute the output based modal amplitudes
    W_meas[i,:] = dot(UT,Z_meas[i,:])

for i in range(nb):
    # read in the results from measured geometries
    rundir = '/mnt/pwfiles/ericdow/random_blade/runs_measured/run'+'%04d/' % (i+1)
    eta_meas[i], pr_meas[i], mdot_meas[i], dbeta_meas[i] = read_utcfd_meas(rundir)

sls_meas = pr_meas - alpha*mdot_meas

pylab.figure()
img = pylab.contourf(W1,W2,eta_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W_meas[:nb,0],W_meas[:nb,1],100,eta_meas,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Efficiency Regression Surface')

pylab.figure()
img = pylab.contourf(W1,W2,sls_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W_meas[:nb,0],W_meas[:nb,1],100,sls_meas,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Speedline Shift Regression Surface')

# compute predicted performance for measured blades
eta_meas_reg = c_C1[0] + c_C1[1]*W_meas[:,0] + c_C1[2]*W_meas[:,1] +\
               c_C1[3]*W_meas[:,0]*W_meas[:,1] + c_C1[4]*W_meas[:,0]**2 + c_C1[5]*W_meas[:,1]**2
sls_meas_reg = c_C2[0] + c_C2[1]*W_meas[:,0] + c_C2[2]*W_meas[:,1] +\
               c_C2[3]*W_meas[:,0]*W_meas[:,1] + c_C2[4]*W_meas[:,0]**2 + c_C2[5]*W_meas[:,1]**2

pylab.figure()
xx = linspace(min(C1),max(C1),2)
pylab.plot(xx,xx,'--')
pylab.plot(eta_meas_reg,eta_meas,'kx',mew=2)
pylab.xlabel('Predicted efficiency')
pylab.ylabel('Actual efficiency')
pylab.grid(True)
pylab.xlim((min(C1),max(C1)))
pylab.ylim((min(C1),max(C1)))

mse_eta_quad = sum((eta_meas_reg-eta_meas)**2)/nb
 
pylab.figure()
xx = linspace(min(C2),max(C2),2)
pylab.plot(xx,xx,'--')
pylab.plot(sls_meas_reg,sls_meas,'kx',mew=2)
pylab.xlabel('Predicted speedline shift')
pylab.ylabel('Actual speedline shift')
pylab.grid(True)
pylab.xlim((min(C2),max(C2)))
pylab.ylim((min(C2),max(C2)))

mse_sls_quad = sum((sls_meas_reg-sls_meas)**2)/nb

# plot histograms of the performance vs the error
eta_samp_reg = c_C1[0] + c_C1[1]*W[:,0] + c_C1[2]*W[:,1] +\
               c_C1[3]*W[:,0]*W[:,1] + c_C1[4]*W[:,0]**2 + c_C1[5]*W[:,1]**2
sls_samp_reg = c_C2[0] + c_C2[1]*W[:,0] + c_C2[2]*W[:,1] +\
               c_C2[3]*W[:,0]*W[:,1] + c_C2[4]*W[:,0]**2 + c_C2[5]*W[:,1]**2

pylab.figure()
pylab.hist(C1-mean(C1),bins=50,facecolor='y')
pylab.hist(eta-eta_samp_reg,bins=50,facecolor='b')
pylab.xlabel('Efficiency')
pylab.ylabel('Frequency')

pylab.figure()
pylab.hist(C2-mean(C2),bins=50,facecolor='y')
pylab.hist(C2-sls_samp_reg,bins=50,facecolor='b')
pylab.xlabel('Speedline shift')
pylab.ylabel('Frequency')

# try a linear regression surface
A = ones((Nsim,3))
Y = copy(C1)
A[:,1] = W[:,0]
A[:,2] = W[:,1]
d_eta = linalg.lstsq(A,Y)[0]

# plot the regression surface
W1 = linspace(W[:,0].min(),W[:,0].max(),100)
W2 = linspace(W[:,1].min(),W[:,1].max(),100)
eta_lin_reg = zeros((100,100))
for i in range(100):
    for j in range(100):
        eta_lin_reg[i,j] = d_eta[0] + d_eta[1]*W1[i] + d_eta[2]*W2[j]
pylab.figure()
img = pylab.contourf(W1,W2,eta_lin_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W[:,0],W[:,1],100,C1,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Efficiency Linear Regression Surface')

A = ones((Nsim,3))
Y = copy(C2)
A[:,1] = W[:,0]
A[:,2] = W[:,1]
d_sls = linalg.lstsq(A,Y)[0]

# plot the regression surface
W1 = linspace(W[:,0].min(),W[:,0].max(),100)
W2 = linspace(W[:,1].min(),W[:,1].max(),100)
sls_lin_reg = zeros((100,100))
for i in range(100):
    for j in range(100):
        sls_lin_reg[i,j] = d_sls[0] + d_sls[1]*W1[i] + d_sls[2]*W2[j]
pylab.figure()
img = pylab.contourf(W1,W2,sls_lin_reg.T,50,cmap=cm.jet)
clim = img.get_clim()
img = pylab.scatter(W[:,0],W[:,1],100,C2,cmap=cm.jet)
img.set_clim(clim)
pylab.xlabel('W1')
pylab.ylabel('W2')
pylab.xlim([W1[0],W1[-1]])
pylab.ylim([W2[0],W2[-1]])
pylab.colorbar()
pylab.title('Speedline Shift Linear Regression Surface')

# plot histograms of the performance vs the error
eta_samp_reg = d_eta[0] + d_eta[1]*W[:,0] + d_eta[2]*W[:,1]
sls_samp_reg = d_sls[0] + d_sls[1]*W[:,0] + d_sls[2]*W[:,1]

pylab.figure()
pylab.hist(C1-mean(C1),bins=50,facecolor='y')
pylab.hist(eta-eta_samp_reg,bins=50,facecolor='b')
pylab.xlabel('Efficiency')
pylab.ylabel('Frequency')

pylab.figure()
pylab.hist(C2-mean(C2),bins=50,facecolor='y')
pylab.hist(C2-sls_samp_reg,bins=50,facecolor='b')
pylab.xlabel('Speedline shift')
pylab.ylabel('Frequency')

# compute predicted performance for measured blades
eta_meas_reg = d_eta[0] + d_eta[1]*W_meas[:,0] + d_eta[2]*W_meas[:,1]
sls_meas_reg = d_sls[0] + d_sls[1]*W_meas[:,0] + d_sls[2]*W_meas[:,1]

pylab.figure()
xx = linspace(min(C1),max(C1),2)
pylab.plot(xx,xx,'--')
pylab.plot(eta_meas_reg,eta_meas,'kx',mew=2)
pylab.xlabel('Predicted efficiency')
pylab.ylabel('Actual efficiency')
pylab.grid(True)
pylab.xlim((min(C1),max(C1)))
pylab.ylim((min(C1),max(C1)))

mse_eta_lin = sum((eta_meas_reg-eta_meas)**2)/nb
 
pylab.figure()
xx = linspace(min(C2),max(C2),2)
pylab.plot(xx,xx,'--')
pylab.plot(sls_meas_reg,sls_meas,'kx',mew=2)
pylab.xlabel('Predicted speedline shift')
pylab.ylabel('Actual speedline shift')
pylab.grid(True)
pylab.xlim((min(C2),max(C2)))
pylab.ylim((min(C2),max(C2)))

mse_sls_lin = sum((sls_meas_reg-sls_meas)**2)/nb
'''

'''    
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
'''

pylab.show()

