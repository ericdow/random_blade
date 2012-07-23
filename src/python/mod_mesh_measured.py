import perturb_surf_measured
import os, shutil

def modify(src,inp,cg_mesh_orig,cg_mesh_mod,rpath,wpath,Z_meas):
    # src          : location of random_blade fortran code
    # inp          : location of CGNS files
    # cg_mesh_orig : name of original CGNS file
    # cg_mesh_mod  : name of modified CGNS file
    # rpath        : name of original blade surface file
    # wpath        : name of modified blade surface file
    # Z_path       : where to store the Z values
    # npca         : number of PCA modes to use in generating geometries

    os.chdir(inp)
    
    # switch to alternate mesh file if original destroyed
    if not os.path.exists(cg_mesh_orig):
        cg_mesh_orig = cg_mesh_orig[:-5]+'1.cgns'
        if not os.path.exists(cg_mesh_orig):
            cg_mesh_orig = cg_mesh_orig[:-5]+'2.cgns'

    # convert CGNS files to HDF format
    os.system('adf2hdf ' + cg_mesh_orig)
    
    # copy CGNS files
    shutil.copy(cg_mesh_orig, cg_mesh_mod)
    
    # write mesh surface out
    os.system(src+'random_blade '+cg_mesh_orig)
    
    # perturb the mesh surface
    perturb_surf_measured.perturb(inp+rpath, inp+wpath, inp+'pca_modes/', Z_meas)
    
    # read in perturbation to CGNS mesh
    os.system(src+'random_blade '+cg_mesh_orig+' '+cg_mesh_mod)
    
    # convert CGNS file back to ADF format
    os.system('hdf2adf ' + cg_mesh_orig)
    os.system('hdf2adf ' + cg_mesh_mod)
