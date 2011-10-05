import perturb_surf
import os, shutil

# file locations
src = '/home/ericdow/code/random_blade/src/'
inp = '/home/ericdow/code/random_blade/input/'

# CGNS mesh
cg_mesh_orig = 'rotor37_fine_inlet.cgns'
cg_mesh_mod = 'rotor37_fine_inlet_mod.cgns'

# CGNS solution
cg_soln_orig = 'design_point_computation_1.cgns'
cg_soln_mod = 'design_point_computation_1_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

os.chdir(inp)

# convert CGNS files to HDF format
os.system('adf2hdf ' + cg_mesh_orig)
# os.system('adf2hdf ' + cg_soln_orig)

# copy CGNS files
shutil.copy(cg_mesh_orig, cg_mesh_mod)
# shutil.copy(cg_soln_orig, cg_soln_mod)

# write mesh surface out
os.system(src+'random_blade '+cg_mesh_orig)
# os.system(src+'random_blade '+cg_soln_orig)

# perturb the mesh surface
perturb_surf.perturb(inp+rpath, inp+wpath)

# read in perturbation to CGNS mesh
os.system(src+'random_blade '+cg_mesh_orig+' '+cg_mesh_mod)
# os.system(src+'random_blade '+cg_soln_orig+' '+cg_soln_mod)

# convert CGNS files back to ADF format
os.system('hdf2adf ' + cg_mesh_mod)
# os.system('hdf2adf ' + cg_soln_mod)
