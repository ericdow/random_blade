import mod_mesh_trans
import os, shutil

ut_src = '/mnt/pwfiles/ericdow/utcfd/bin/'
src = '/mnt/pwfiles/ericdow/random_blade/src/'
inp = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/'

cg_mesh_orig = 'utcfd_out.cgns'
cg_mesh_mod  = 'utcfd_out_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

# number of iterations to run each mesh
niter = 2000
# number of cores to run on
np = 1
# mode to scale
mode = 2
# scale of the modal perturbation
scale = -2.0

rundir = '/mnt/pwfiles/ericdow/random_blade/trans_test/mode'+str(mode)+'/scale'+str(scale)+'/'
if not(os.path.exists(rundir)):
    os.mkdir(rundir)

# generate modified mesh
mod_mesh_trans.modify(src,inp,cg_mesh_orig,cg_mesh_mod,rpath,wpath,scale,mode)
shutil.move(inp+cg_mesh_mod,rundir+'utcfd_out_mod.cgns')

# copy blade_surf_mod.dat file
shutil.copy(inp+'blade_surf_mod.dat',rundir)

# copy casefile to run directory
# shutil.copy(inp+'casefile.inp',rundir)
shutil.copy(inp+'cart2cyl.inp',rundir)

# convert mesh to cylindrical 
os.chdir(rundir)
os.system(ut_src+'preut -f utcfd_out_mod.cgns -o utcfd_in -b cart2cyl.inp')

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

