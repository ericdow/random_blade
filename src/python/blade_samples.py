import perturb_surf_gaussian
import os, shutil

ut_src = '/mnt/pwfiles/ericdow/utcfd/bin/'
src = '/mnt/pwfiles/ericdow/random_blade/src/'
inp = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/'
runs = '/mnt/pwfiles/ericdow/random_blade/runs/'

cg_mesh_orig = 'utcfd_out.cgns'
cg_mesh_mod  = 'utcfd_out_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

# number of KL modes to use in generating geometries
nkl = 50

# perturb the mesh surface
perturb_surf_gaussian.perturb(inp+rpath, 'foo.dat', inp+'kl_modes/', 'foo.dat', 'foo.dat', nkl)
