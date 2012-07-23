import perturb_surf_pca
import os, shutil
from filelock import FileLock

def modify(src,inp,rundir,cg_mesh_orig,cg_mesh_mod,rpath,wpath,Z_path,npca):
    # src          : location of random_blade fortran code
    # inp          : location of CGNS files
    # rundir       : current running directory
    # cg_mesh_orig : name of original CGNS file
    # cg_mesh_mod  : name of modified CGNS file
    # rpath        : name of original blade surface file
    # wpath        : name of modified blade surface file
    # Z_path       : where to store the Z values
    # npca         : number of PCA modes to use in generating geometries

    with FileLock(inp+cg_mesh_orig,timeout=20):
        # copy CGNS to rundir
        shutil.copy(inp+cg_mesh_orig,rundir)

    os.chdir(rundir)

    # convert CGNS files to HDF format
    os.system('adf2hdf ' + cg_mesh_orig)
    
    # copy CGNS file
    shutil.copy(cg_mesh_orig, cg_mesh_mod)
    
    # write mesh surface out
    os.system(src+'random_blade '+cg_mesh_orig)
    
    # perturb the mesh surface
    perturb_surf_pca.perturb(inp+rpath, rundir+wpath, inp+'pca_modes_gp/', Z_path, npca)
    
    # read in perturbation to CGNS mesh
    os.system(src+'random_blade '+cg_mesh_orig+' '+cg_mesh_mod)
    
    # convert CGNS file back to ADF format
    os.system('hdf2adf ' + cg_mesh_mod)

    # remove the original CGNS file
    os.remove(cg_mesh_orig)
