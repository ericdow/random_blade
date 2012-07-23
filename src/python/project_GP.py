import pylab, os, math, string
from numpy import *
from numpy.linalg import *
from matplotlib import cm
from scipy.interpolate import splprep, splev, sproot, splrep

def stretch_func(t,slope):
    # slope: slop at t = 0,1
    A = ones((3,3))
    A[1,:2] = 0.0
    A[2,0] = 3.0
    A[2,1] = 2.0
    b = array((1.0,slope,slope))
    abc = linalg.solve(A,b)

    return abc[0]*t**3 + abc[1]*t**2 + abc[2]*t

def project(s,t,sABCD,fp,x,y,z,slope):

    # interpolate Gaussian random field to mesh

    # inputs
    # s,t   : coordinates of the grid on which GP is simulated
    # sABCD : s-coordinate of the map control points on GP mesh
    # fp    : values of the GP on the grid
    # xyz   : coordinates of blade to transform to (from mesh)

    cdim = x.shape[0]
    sdim = x.shape[1]

    # normalize
    s /= s[-1]
    t /= t[-1] 

    '''
    for i in range(len(s)):
        if (s[i] < sABCD[0]) or (s[i] > sABCD[3]):
            fp[i,:] = 0.0
        if (s[i] >= sABCD[0]) and (s[i] < sABCD[1]):
            fp[i,:] = 1.0
        if (s[i] >= sABCD[1]) and (s[i] < sABCD[2]):
            fp[i,:] = 0.0
        if (s[i] >= sABCD[2]) and (s[i] < sABCD[3]):
            fp[i,:] = 1.0
    '''

    fp_mesh = zeros((cdim,sdim))

    # create a spline representation of each section of the random field
    # fp_interp[3*i]   - t for ith section
    # fp_interp[3*i+1] - c for ith section
    # fp_interp[3*i+2] - k for ith section
    ns = len(s)
    nt = len(t)
    for i in range(nt):
        tck_tmp = splrep(s,fp[:,i])
        if i == 0:
            fp_interp = tck_tmp
        else:
            fp_interp = fp_interp+tck_tmp

    tck_LE, t_LE = splprep([x[(cdim+1)/2,:],y[(cdim+1)/2,:],z[(cdim+1)/2,:]],s=0)
    tck_TE, t_TE = splprep([x[0,:],y[0,:],z[0,:]],s=0)

    # create a spline representation of each spanwise section of the mesh
    s_mesh = zeros((sdim,cdim))
    s_LE = zeros((sdim))
    s_A = zeros((sdim))
    s_B = zeros((sdim))
    s_C = zeros((sdim))
    s_D = zeros((sdim))
    for isec in range(sdim):
        tck_tmp,s_mesh[isec,:] = splprep([x[:,isec],y[:,isec],z[:,isec]],s=0,per=1)
        if isec == 0:
            tck_s = tck_tmp
        else:
            tck_s = tck_s+tck_tmp

        s_LE[isec] = s_mesh[isec,(cdim+1)/2]

        '''
        # TODO: for now, points A,B,C,D are a fixed percent of the blade circumference
        perc = 0.05;
        s_A[isec] = perc
        s_B[isec] = s_LE[isec]-perc
        s_C[isec] = s_LE[isec]+perc
        s_D[isec] = 1.0 - perc
        '''
        s_A[isec] = 0.00208224 + (0.0016914 - 0.00208224)*t_TE[isec]
        s_B[isec] = (s_LE[isec]-0.00352707) + \
                   ((s_LE[isec]-0.00145799)-(s_LE[isec]-0.00352707))*t_LE[isec]
        s_C[isec] = (s_LE[isec]+0.00248137) + \
                   ((s_LE[isec]+0.00191383)-(s_LE[isec]+0.00248137))*t_LE[isec]
        s_D[isec] = 0.997917 + (0.998497 - 0.997917)*t_TE[isec]
                
    # stretch things towards the hub and shroud
    tstretch = stretch_func(t,slope)

    # perform the interpolation
    for isec in range(sdim):
        for i in range(cdim):

            # interpolate region from 0 to s_A
            if (s_mesh[isec,i] <= s_A[isec]):
                ss = (s_mesh[isec,i]+1.0-s_D[isec])/(s_A[isec]+1.0-s_D[isec])
                # find the t coordinate of this point
                xyz = zeros((sdim,3))
                for ii in range(sdim):
                    sAD = s_A[ii]+1.0-s_D[ii]
                    if ss >= (1.0-s_D[ii])/sAD:
                        sev = ss*sAD-(1.0-s_D[ii])
                    else:
                        sev = s_D[ii] + ss*sAD
                    # spline representation of this spanwise section
                    tck = (tck_s[3*ii],tck_s[3*ii+1],tck_s[3*ii+2])
                    xyz[ii,0],xyz[ii,1],xyz[ii,2] = splev(sev,tck)
                foo,tt = splprep([xyz[:,0],xyz[:,1],xyz[:,2]],s=0,per=0)
                # s coordinate on the GP mesh
                sAD = sABCD[0]+1.0-sABCD[3]
                if ss >= (1.0-sABCD[3])/sAD:
                    sGP = ss*sAD-(1.0-sABCD[3])
                else:
                    sGP = sABCD[3] + ss*sAD
                fp_t = zeros(nt)
                for ii in range(nt):
                    # spline representation of this spanwise section
                    tck = (fp_interp[3*ii],fp_interp[3*ii+1],fp_interp[3*ii+2])
                    fp_t[ii] = splev(sGP,tck)
                tck = splrep(tstretch,fp_t)
                fp_mesh[i,isec] = splev(tt[isec],tck)

            # interpolate region from s_D to 0
            if (s_mesh[isec,i] >= s_D[isec]):
                ss = (s_mesh[isec,i]-s_D[isec])/(s_A[isec]+1.0-s_D[isec])
                # find the t coordinate of this point
                xyz = zeros((sdim,3))
                for ii in range(sdim):
                    sAD = s_A[ii]+1.0-s_D[ii]
                    if ss >= (1.0-s_D[ii])/sAD:
                        sev = ss*sAD-(1.0-s_D[ii])
                    else:
                        sev = s_D[ii] + ss*sAD
                    # spline representation of this spanwise section
                    tck = (tck_s[3*ii],tck_s[3*ii+1],tck_s[3*ii+2])
                    xyz[ii,0],xyz[ii,1],xyz[ii,2] = splev(sev,tck)
                foo,tt = splprep([xyz[:,0],xyz[:,1],xyz[:,2]],s=0,per=0)
                # s coordinate on the GP mesh
                sAD = sABCD[0]+1.0-sABCD[3]
                if ss >= (1.0-sABCD[3])/sAD:
                    sGP = ss*sAD-(1.0-sABCD[3])
                else:
                    sGP = sABCD[3] + ss*sAD
                fp_t = zeros(nt)
                for ii in range(nt):
                    # spline representation of this spanwise section
                    tck = (fp_interp[3*ii],fp_interp[3*ii+1],fp_interp[3*ii+2])
                    fp_t[ii] = splev(sGP,tck)
                tck = splrep(tstretch,fp_t)
                fp_mesh[i,isec] = splev(tt[isec],tck)

            # interpolate region from s_A to s_B
            if (s_mesh[isec,i] > s_A[isec]) and (s_mesh[isec,i] <= s_B[isec]):
                ss = (s_mesh[isec,i]-s_A[isec])/(s_B[isec]-s_A[isec])
                # find the t coordinate of this point
                xyz = zeros((sdim,3))
                for ii in range(sdim):
                    sAB = s_B[ii]-s_A[ii]
                    # spline representation of this spanwise section
                    tck = (tck_s[3*ii],tck_s[3*ii+1],tck_s[3*ii+2])
                    xyz[ii,0],xyz[ii,1],xyz[ii,2] = splev(sev,tck)
                foo,tt = splprep([xyz[:,0],xyz[:,1],xyz[:,2]],s=0,per=0)
                # s coordinate on the GP mesh
                sAB = sABCD[1] - sABCD[0]
                sGP = sABCD[0] + ss*sAB
                fp_t = zeros(nt)
                for ii in range(nt):
                    # spline representation of this spanwise section
                    tck = (fp_interp[3*ii],fp_interp[3*ii+1],fp_interp[3*ii+2])
                    fp_t[ii] = splev(sGP,tck)
                tck = splrep(tstretch,fp_t)
                fp_mesh[i,isec] = splev(tt[isec],tck)

            # interpolate region from s_B to s_C
            if (s_mesh[isec,i] >= s_B[isec]) and (s_mesh[isec,i] < s_C[isec]):
                ss = (s_mesh[isec,i]-s_B[isec])/(s_C[isec]-s_B[isec])
                # find the t coordinate of this point
                xyz = zeros((sdim,3))
                for ii in range(sdim):
                    sBC = s_C[ii]-s_B[ii]
                    # spline representation of this spanwise section
                    tck = (tck_s[3*ii],tck_s[3*ii+1],tck_s[3*ii+2])
                    xyz[ii,0],xyz[ii,1],xyz[ii,2] = splev(sev,tck)
                foo,tt = splprep([xyz[:,0],xyz[:,1],xyz[:,2]],s=0,per=0)
                # s coordinate on the GP mesh
                sBC = sABCD[2] - sABCD[1]
                sGP = sABCD[1] + ss*sBC
                fp_t = zeros(nt)
                for ii in range(nt):
                    # spline representation of this spanwise section
                    tck = (fp_interp[3*ii],fp_interp[3*ii+1],fp_interp[3*ii+2])
                    fp_t[ii] = splev(sGP,tck)
                tck = splrep(tstretch,fp_t)
                fp_mesh[i,isec] = splev(tt[isec],tck)

            # interpolate region from s_C to s_D
            if (s_mesh[isec,i] >= s_C[isec]) and (s_mesh[isec,i] < s_D[isec]):
                ss = (s_mesh[isec,i]-s_C[isec])/(s_D[isec]-s_C[isec])
                # find the t coordinate of this point
                xyz = zeros((sdim,3))
                for ii in range(sdim):
                    sCD = s_D[ii]-s_C[ii]
                    # spline representation of this spanwise section
                    tck = (tck_s[3*ii],tck_s[3*ii+1],tck_s[3*ii+2])
                    xyz[ii,0],xyz[ii,1],xyz[ii,2] = splev(sev,tck)
                foo,tt = splprep([xyz[:,0],xyz[:,1],xyz[:,2]],s=0,per=0)
                # s coordinate on the GP mesh
                sCD = sABCD[3] - sABCD[2]
                sGP = sABCD[2] + ss*sCD
                fp_t = zeros(nt)
                for ii in range(nt):
                    # spline representation of this spanwise section
                    tck = (fp_interp[3*ii],fp_interp[3*ii+1],fp_interp[3*ii+2])
                    fp_t[ii] = splev(sGP,tck)
                tck = splrep(tstretch,fp_t)
                fp_mesh[i,isec] = splev(tt[isec],tck)

    # make sure things match periodically
    fp_mesh[0,:] = fp_mesh[-1,:]

    return fp_mesh

