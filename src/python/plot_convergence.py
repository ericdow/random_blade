import pylab, os, math, string
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import read_blade, write_tecplot

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
                nsteps = int(lines[-4].strip().split()[0])
            lines = file(rundir+str(rd)+'/Z.dat').readlines()
            npca = len(lines)
        i += 1

    nruns = int(sum(iran))
    resid = zeros((nruns,nsteps))

    # find which lines to read in
    ir = 0
    js = -1
    for rd in os.listdir(rundir):
        if iran[ir] and js == -1:
            lines = file(rundir+str(rd)+'/utcfd_out.monitoring.dat').readlines()
            j = 0
            for line in lines:
                if line[:15] == ' Job is started':
                    js = j
                j += 1
        ir += 1

    ir = 0
    i = 0
    for rd in os.listdir(rundir):
        if iran[ir]:
            lines = file(rundir+str(rd)+'/utcfd_out.monitoring.dat').readlines()
            # check for negative volumes
            check = filter(lambda x: x[:25] == ' Warning: Negative volume',lines)
            if check:
                print rd
            else:
                lines = lines[js+2:-3]
                # remove check point lines
                lines = filter(lambda x: x[0] != 'C',lines)
                j = 0
                for line in lines:
                    tmp = line.strip().split()
                    resid[i,j] = tmp[2]
                    j += 1
                i += 1
        ir += 1

    return nsteps, resid

# read in the results
rundir = '/mnt/pwfiles/ericdow/random_blade/runs/'
nsteps, resid = read_utcfd(rundir)

pylab.plot(arange(nsteps),resid[:100,:].T)
pylab.show()


