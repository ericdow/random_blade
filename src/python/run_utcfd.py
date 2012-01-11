import perturb_surf
import mod_mesh
import os, shutil

ut_src = '/mnt/pwfiles/ericdow/utcfd/bin/'
src = '/home/ericdow/code/random_blade/src/'
inp = '/home/ericdow/code/random_blade/input/'
runs = '/home/ericdow/code/random_blade/runs/'

cg_mesh_orig = 'utcfd_out.cgns'
cg_mesh_mod  = 'utcfd_out_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

# number of MC samples to run
n_mc = 1
# number of iterations to run each mesh
niter = 1000
# number of cores to run on
np = 1

for i in range(n_mc):
    # generate modified mesh
    mod_mesh.modify(src,inp,cg_mesh_orig,cg_mesh_mod,rpath,wpath)

    # create directory to run utcfd in
    rundir = runs+'run'+'%04d' % i)
    os.mkdir(runsdir)
    shutil.move(inp+cg_mesh_mod,rundir+'utcfd_in.cgns')

    # run preut
    os.chdir(rundir)
    os.system(ut_src+'preut -f utcfd_in.cgns -b casefile.inp')

    # copy correct 'utcfd.bc' to run directory
    shutil.copy(runs+'utcfd.bc',rundir)

    # run utcfd
    os.system('mpiexec -np '+str(np)+' 'ut_src+'utcfd.exe'

    # process results
