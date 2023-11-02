# superficial-White-Matter
Generates surfaces at various white matter depths (default 1, 2, and 3%)
![example](https://github.com/jordandekraker/superficial-white-matter/blob/main/scrnshot.png)
Red is the original wm surface, orange-yellow are depths 1mm, 2mm, and 3mm

This is done by first computing a Laplace field over white matter (cortex to subcortex+ventricles), and then shifting an exiting white matter surface along that gradient. Stopping conditions are set by geodesic distance travelled.

## Installation
```
git clone https://github.com/jordandekraker/superficial-white-matter.git
pip install superficial-white-matter/
```

## Usage (example)
```
python sWM/laplace_solver.py \
  EXAMPLE/sub-01_aparc+aseg_space-nativepro.nii.gz \
  OUT/sub-01_space-nativepro_laplace-wm.nii.gz \
python sWM/surface_generator.py \
  EXAMPLE/sub-01_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii \
  OUT/sub-01_space-nativepro_laplace-wm.nii.gz \
  OUT/sub-01_hemi-R_space-nativepro_surf-fsLR-32k_label-sWF_depth-
```
