#/bin/bash

# Solve a Laplace field
python laplace_solver.py aparc+aseg_space-nativepro.nii.gz wm-laplace.nii.gz

# Shift a given surface along the Laplace field
for depths in 0.05 0.1 0.2 0.3; do
    python laplace_surf_interp.py sub-01_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii wm-laplace.nii.gz wm-equipotentials/sub-01_hemi-L_space-nativepro_surf-fsLR-32k_SFWdepth-${depths}" $depths
done
