#/bin/bash

DERIV='/data/mica3/BIDS_MICs/derivatives'
SUB='HC002'
SES='01'

# import some nifti/gifti data
mkdir test
mri_convert $DERIV/fastsurfer/sub-${SUB}_ses-${SES}/mri/aparc+aseg.mgz test/aparc+aseg.nii.gz
tform='' #$DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/xfm/sub-${SUB}_ses-${SES}_from-fsnative_to_nativepro_T1w_0GenericAffine.mat
antsApplyTransforms -i test/aparc+aseg.nii.gz -r $DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/anat/sub-${SUB}_ses-${SES}_space-nativepro_T1w.nii.gz -o test/aparc+aseg_space-nativepro.nii.gz -n MultiLabel
cp $DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/surf/sub-${SUB}_ses-${SES}_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii test/wm.surf.gii

# run
python sfw/laplace_solver.py test/aparc+aseg_space-nativepro.nii.gz test/wm-laplace.nii.gz
python sfw/surface_generator.py test/wm.surf.gii test/wm-laplace.nii.gz test/depth
