import perturb_surf_gaussian
import mod_mesh_kl
import os, shutil

src = '/mnt/pwfiles/ericdow/random_blade/src/'
inp = '/mnt/pwfiles/ericdow/random_blade/input/rotor37_coarse/'
runs = '/mt/pwfiles/ericdow/random_blade/runs/'

cg_mesh_orig = 'utcfd_out.cgns'
cg_mesh_mod  = 'utcfd_out_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

perturb_surf_gaussian.perturb(rpath, wpath)

'''
# 33 total PCA modes
npca = 33
perturb_surf_pca.perturb(rpath, wpath, inp+'pca/','Z.dat',npca)
'''
