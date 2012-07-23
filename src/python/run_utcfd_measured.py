import mod_mesh_measured
import os, shutil
from numpy import *

ut_src = '/mnt/pwfiles/ericdow/utcfd/bin/'
src = '/mnt/pwfiles/ericdow/random_blade/src/'
inp = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/'
runs = '/mnt/pwfiles/ericdow/random_blade/runs_measured/'

cg_mesh_orig = 'utcfd_out.cgns'
cg_mesh_mod  = 'utcfd_out_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

# number of blades to run
n_blades = 34
# number of iterations to run each mesh
niter = 3000
# number of cores to run on
np = 1

Z_meas = zeros((n_blades,n_blades))

for iblade in range(n_blades):
    # create directory to run utcfd in
    rundir = runs+'run'+'%04d/' % (iblade+1)
    os.mkdir(rundir)

    # read in the Z values for each blade
    lines = file('/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/pca_modes/Z_blade'+str(iblade+1)+'.dat').readlines()
    Z_meas[iblade,:] = array([line.strip().split() for line in lines], float).squeeze()

    # generate modified mesh
    meas_path = 'measured_blades/blade'+str(iblade+1)+'.dat'
    mod_mesh_measured.modify(src,inp,cg_mesh_orig,cg_mesh_mod,rpath,wpath,Z_meas[iblade,:])
    shutil.move(inp+cg_mesh_mod,rundir+'utcfd_out_mod.cgns')

    # copy blade_surf_mod.dat file
    shutil.copy(inp+'blade_surf_mod.dat',rundir)

    # copy casefile to run directory
    # shutil.copy(inp+'casefile.inp',rundir)
    shutil.copy(inp+'cart2cyl.inp',rundir)

    # convert mesh to cylindrical 
    os.chdir(rundir)
    os.system(ut_src+'preut -f utcfd_out_mod.cgns -o utcfd_in -b cart2cyl.inp')

    # remove the modified mesh file
    os.remove('utcfd_out_mod.cgns')

    # run preut (file must be in adf format to run preut)
    # os.system(ut_src+'preut -f utcfd_in.cgns -b casefile.inp')

    # copy correct 'utcfd.bc' to run directory
    shutil.copy(inp+'utcfd.bc',rundir)

    # set number of iterations in 'utcfd.input'
    lines = file('utcfd.input').readlines()
    il = 0
    for line in lines:
        line = line.strip().split()
        if line:
            if line[0] == 'nstep':
                lines[il] = lines[il].replace(line[-1],str(niter))
        il += 1
    
    f = open('utcfd.input','w')
    f.writelines(lines)
    f.close()

    # run utcfd
    os.system('mpiexec -np '+str(np)+' '+ut_src+'utcfd.exe > out')

