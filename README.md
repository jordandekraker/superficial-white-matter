# superficial-White-Matter
Generates surfaces at various white matter depths (default 1,2, and 3%)
![example](https://github.com/jordandekraker/superficial-white-matter/blob/main/scrnshot.png)
Red isthe original wm surface, orange-yellow are 1-3% depths (default)

This is done by first computing a Laplace field over white matter (between cortex to subcortex), and then shifting an exiting white matter surface along that gradient. Stopping conditions are set by thresholding the Laplace field.

## Installation
'''
git clone https://github.com/jordandekraker/superficial-white-matter.git
cd superficial-white-matter
pip install .
'''

## Usage (example)
'''
python sWM/laplace_solver.py EXAMPLE/sub-01_aparc+aseg_space-nativepro.nii.gz OUT/sub-01_wm-laplace.nii.gz
python sWM/surface_generator.py EXAMPLE/sub-01_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii OUT/sub-01_wm-laplace.nii.gz OUT/sub-01_hemi-R_space-nativepro_surf-fsLR-32k_label-sWF_depth-
'''
