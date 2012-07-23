import pylab, os, math, string, shutil
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# read in the results
rundir = '/mnt/pwfiles/ericdow/random_blade/runs_choke_mdot/'

# remove empty/unrun directories
iran = zeros((len(os.listdir(rundir))))
i = 0
for rd in os.listdir(rundir):
    # delete if output directory doesn't exist
    if str(rd[0:3]) == 'run':
        if not(os.path.exists(rundir+str(rd)+'/output')):
            print 'deleting '+rundir+str(rd)
            shutil.rmtree(rundir+str(rd))
        # delete if output directory is empty
        if os.path.exists(rundir+str(rd)+'/output'):
            if os.listdir(rundir+str(rd)+'/output') == []:
                print 'deleting '+rundir+str(rd)
                shutil.rmtree(rundir+str(rd))

# remove utcfd_out_mod CGNS file
for rd in os.listdir(rundir):
    if os.path.exists(rundir+str(rd)+'/utcfd_out_mod.cgns'):
        os.remove(rundir+str(rd)+'/utcfd_out_mod.cgns')

