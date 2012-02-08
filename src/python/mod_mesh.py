import perturb_surf
import os, shutil

def modify(src,inp,cg_mesh_orig,cg_mesh_mod,rpath,wpath):
    # src          : location of random_blade fortran code
    # inp          : location of CGNS files
    # cg_mesh_orig : name of original CGNS file
    # cg_mesh_mod  : name of modified CGNS file
    # rpath        : name of original blade surface file
    # wpath        : name of modified blade surface file

    os.chdir(inp)
    
    # convert CGNS files to HDF format
    os.system('adf2hdf ' + cg_mesh_orig)
    
    # copy CGNS files
    shutil.copy(cg_mesh_orig, cg_mesh_mod)
    
    # write mesh surface out
    os.system(src+'random_blade '+cg_mesh_orig)
    
    # perturb the mesh surface
    perturb_surf.perturb(inp+rpath, inp+wpath)
    
    # read in perturbation to CGNS mesh
    os.system(src+'random_blade '+cg_mesh_orig+' '+cg_mesh_mod)
    
    # convert CGNS file back to ADF format
    os.system('hdf2adf ' + cg_mesh_orig)
    os.system('hdf2adf ' + cg_mesh_mod)
