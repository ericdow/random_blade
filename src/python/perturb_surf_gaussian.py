import read_blade, simulate, write_tecplot, project_GP
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def perturb(rpath, wpath, kl_path, Z_path, dth_path, nkl):
    
    # inputs
    # Z_path   : where to write out the sample normal 
    # dth_path : where to write out the random rotation
    # nkl      : where to truncate the K-L expansion

    cdim, sdim, x, y, z = read_blade.read_coords(rpath)
    
    # translate blade random amount (rotate about x axis)
    random.seed()
    dth = random.randn(1)*pi/180*0.1/3.0
    for i in range(cdim):
        for j in range(sdim):
            r = sqrt(y[i,j]**2 + z[i,j]**2)
            th = math.atan2(z[i,j],y[i,j])
            y[i,j] = r*cos(th+dth)
            z[i,j] = r*sin(th+dth)

    # write out dth
    f = open(dth_path,'w')
    f.write('%e\n' % dth)
    f.close()
    
    # generate and write out the Z for this mesh
    random.seed() 
    Z = random.randn(nkl)    # normal distribution vector
    f = open(Z_path,'w')
    for i in range(nkl):
        f.write('%e\n' % Z[i])
    f.close()

    # read in PCA modes and singular values 
    lines = file(kl_path+'S.dat').readlines()
    S = array([line.strip().split()[0] for line in lines], float).T
    n_modes = len(S)
    fp = zeros((cdim,sdim))
    for i in range(nkl):
        cdim,sdim,V = read_blade.read_mode(kl_path+'V'+str(i+1)+'.dat')
        fp += Z[i]*sqrt(S[i])*V
     
    np = read_blade.calcNormals(x,y,z)
    
    xp = x + np[:,:,0]*fp
    yp = y + np[:,:,1]*fp
    zp = z + np[:,:,2]*fp
     
    # write out the blade surface
    f = open(wpath,'w')
    f.write('CDIM:       %d\n' % cdim)
    f.write('SDIM:       %d\n' % sdim)
    for i in arange(cdim):
        for j in arange(sdim):
            f.write('%20.8E' * 3 % (xp[i,j],yp[i,j],zp[i,j]))
            f.write('\n')
    
    f.close()

    # write out to tecplot format
    write_tecplot.write_blade_surf(xp,yp,zp,fp,'blade.dat')

