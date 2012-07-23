import read_blade
import simulate
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def perturb(rpath, wpath, pca_path, scale, mode):

    # inputs
    # Z_path : where to write out the sample normal 
    # npca   : where to truncate the PCA expansion

    cdim, sdim, x, y, z = read_blade.read_coords(rpath)

    # perturb using only one of the transformed modes
    lines = file(pca_path+'S.dat').readlines()
    S = array([line.strip().split()[0] for line in lines], float).T
    fp = zeros((cdim,sdim))
    cdim,sdim,V = read_blade.read_mode(pca_path+'V'+str(mode)+'.dat')
    fp += scale*S[mode-1]*V
     
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

    # test to see how they look
    '''
    i = 3
    cdim,sdim,V = read_blade.read_mode(pca_path+'V'+str(i+1)+'.dat')
    s, t = read_blade.xyz2st(x,y,z)
    print s[0],s[1]
    pylab.contourf(s,t,V.T,50)
    pylab.colorbar()
    pylab.show()
    '''
     
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
       



