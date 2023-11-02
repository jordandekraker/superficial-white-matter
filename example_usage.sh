#/bin/bash

# Assuming you already have micapipe_v0.2.0 and free/fastSurfer

DERIV='/data/mica3/BIDS_MICs/derivatives'
SUB='HC002'
SES='01'

# format some nifti/gifti data
mkdir -p $DERIV/sWM/sub-${SUB}/ses-${SES}
mri_convert \
  $DERIV/fastsurfer/sub-${SUB}_ses-${SES}/mri/aparc+aseg.mgz \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-fsnative_aparc+aseg.nii.gz
tform='' # $DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/xfm/sub-${SUB}_ses-${SES}_from-fsnative_to_nativepro_T1w_0GenericAffine.mat
antsApplyTransforms \
  -i $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-fsnative_aparc+aseg.nii.gz \
  -r $DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/anat/sub-${SUB}_ses-${SES}_space-nativepro_T1w.nii.gz \
  -o $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-nativepro_aparc+aseg.nii.gz \
  -n MultiLabel
#  -t $tform
# even if the tform is not needed, this resampling is necessary to standardize nifti headers

# run
python sWM/laplace_solver.py \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-nativepro_aparc+aseg.nii.gz \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-nativepro_laplace-wm.nii.gz 
python sWM/surface_generator.py \
  $DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/surf/sub-${SUB}_ses-${SES}_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-nativepro_laplace-wm.nii.gz \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_hemi-R_space-nativepro_label-sWF_depth- 
python sWM/surface_generator.py \
  $DERIV/micapipe_v0.2.0/sub-${SUB}/ses-${SES}/surf/sub-${SUB}_ses-${SES}_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_space-nativepro_laplace-wm.nii.gz \
  $DERIV/sWM/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_hemi-L_space-nativepro_label-sWF_depth- 
