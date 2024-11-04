import nibabel as nib
import numpy as np
from scipy.stats import mode
import subprocess
import os
from .utilities import show_warning
from .sWM import laplace_solver, surface_generator
import ants
from ants.ops import resample_image_to_target

def load_gifti_data(filepath):
    data = nib.load(filepath)
    return data.darrays[0].data


def calcdist(surf1, surf2):

    euclidianDistanceSquared = (surf1 - surf2) ** 2
    euclidianDistanceSummed = np.sum(euclidianDistanceSquared, axis=1)
    return np.sqrt(euclidianDistanceSummed)


def computegrad(data, dists):
    data = np.ediff1d(data)
    data[dists == 0] = 0
    dists[dists == 0] = 1
    return np.divide(data, dists)

def maptosurface(mri_map,surf_fs,map_on_surf,hemi,label_data,out_map,idBIDS,workbench_path):
          # -----------------------------------------------
  # Volume to surface mapping
  # Function that maps a MRI volume to a surfaces
  # Using work bench commands and multiple surfaces:
  # fsnative, fsaverage5, fsLR-5k and fsLR-32k
  # -----------------------------------------------
  # Input variables
  #mri_map=$1                      # MRI map from where data will be mapped
  #surf_fs=$2                      # Surface to map the MRI (MUST be surf-fsnative, same space ast mri_map)
  #map_on_surf=${3}                # Outname of the data mapped on the surface
  #H=$4                            # Hemisphere {L, R}
  #label_data=$5                   # label of the map (e.g. FA, ADC, flair, T2star, MTR)
  #out_map=$6
  # Map to highest resolution surface (fsnative: more vertices)
    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "volume-to-surface-mapping",
            mri_map,
            surf_fs,
            map_on_surf,
            "-trilinear",
        ]
    )


def compute_blurring(
    input_dir,
    surf_dir,
    bids_id,
    hemi,
    output_file,
    feat,
    workbench_path,
    resol,
    fwhm,
    surf_file,
    output_file_final,
    tmp_dir,
):

    base_path = input_dir
    for _ in range(4):  # Adjust the range to navigate up the desired number of levels
        base_path, _ = os.path.split(base_path)


    micapipe_path = os.path.split(input_dir)[0]

    freesurfer_path = os.path.join(
        base_path, "freesurfer", bids_id, "mri", "aparc+aseg.nii.gz"
    )
    temp_parc_path = os.path.join(tmp_dir, f"{hemi}_surf-fsnative_label-temp.nii.gz")

    output_path = os.path.join(tmp_dir, f"laplace{hemi}.nii.gz")
    # img = ants.image_read(freesurfer_path)
    # imgfixed = ants.image_read(f"{input_dir}/{bids_id}_space-nativepro_map-T1map.nii.gz")
    # resample_image_to_target(imgfixed, img, interp_type='multiLabel', verbose=True).to_filename(temp_parc_path)
    if not os.path.exists(freesurfer_path):
        return 'Freesurfer path does not exist'
    if not os.path.exists(temp_parc_path):
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-volume-resample",
                freesurfer_path,
                f"{input_dir}/{bids_id}_space-nativepro_map-T1map.nii.gz",
                "ENCLOSING_VOXEL",
                temp_parc_path,
                # "-affine",
                # os.path.join(micapipe_path,"xfm", f"{bids_id}_from-fsnative_to_nativepro_T1w_0GenericAffine.mat")
            ]
        )




        if not os.path.exists(os.path.join(tmp_dir, "swm")):
            os.mkdir(os.path.join(tmp_dir, "swm"))
        laplace_solver.solve_laplace(temp_parc_path, output_path)
        surface_generator.shift_surface(
            f"{surf_dir}/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii",
            output_path,
            f"{tmp_dir}//swm//{hemi}_sfwm-",
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        )

    wmBoundaryDataArr = load_gifti_data(
        f"{input_dir}/{bids_id}_hemi-{hemi}_surf-fsnative_label-white_{feat}.func.gii"
    )
    wmBoundarySurfaceArr = load_gifti_data(
        f"{surf_dir}/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii"
    )

    modeofboundary = mode(wmBoundaryDataArr, keepdims=True)

    midthicknessDataArr = load_gifti_data(
        f"{input_dir}/{bids_id}_hemi-{hemi}_surf-fsnative_label-midthickness_{feat}.func.gii"
    )
    midthicknessSurfaceArr = load_gifti_data(
        f"{surf_dir}/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-midthickness.surf.gii"
    )

    surfarr = [
        [midthicknessDataArr, midthicknessSurfaceArr],
        [wmBoundaryDataArr, wmBoundarySurfaceArr],
    ]
    if feat.lower() != "adc" or feat.lower() != "fa":
        volumemap = f"{input_dir}/{bids_id}_space-nativepro_map-{feat}.nii.gz"
    else:
        volumemap = f"{input_dir}/{bids_id}_space-nativepro_model-DTI_map-{feat.upper()}.nii.gz"
    for surf in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-volume-to-surface-mapping",
                volumemap,
                f"{tmp_dir}//swm//{hemi}_sfwm-{surf}.0mm.surf.gii",
                f"{tmp_dir}//swm//{hemi}_{feat}_{resol}_{fwhm}sfwm-{surf}.0mm-metric.func.gii",
                "-trilinear",
            ]
        )
        surfarr.append([
            load_gifti_data(f"{tmp_dir}//swm//{hemi}_{feat}_{resol}_{fwhm}sfwm-{surf}.0mm-metric.func.gii"),
            load_gifti_data(f"{tmp_dir}//swm//{hemi}_sfwm-{surf}.0mm.surf.gii")
            ])
    


    distances = np.zeros(shape=(len(midthicknessDataArr), len(surfarr) - 1))
    dataArr = np.zeros(shape=(len(midthicknessDataArr), len(surfarr)))
    for e, ds in enumerate(surfarr):
        data, surf = ds
        dataArr[:, e] = np.divide(data, modeofboundary.mode[0])
        if e == len(surfarr) - 1:
            break
        nextdata, nextsurt = surfarr[e + 1]
        distance = calcdist(surf, nextsurt)

        distances[:, e] = distance

    blurring = np.zeros(
        shape=(len(midthicknessDataArr), len(surfarr) - 1), dtype=np.float32
    )
    for i in range(len(dataArr) - 1):
        gradient = computegrad(dataArr[i], distances[i])
        gradient = np.nan_to_num(gradient)
        if gradient[0] == 0:
            gradient = np.zeros_like(gradient)
        blurring[i] = gradient

    # for e, i in enumerate(["midthickness-0mm", "0mm-1mm", "1mm-2mm", "2mm-3mm"]):
    data_array = nib.gifti.gifti.GiftiDataArray(
        data=blurring,
        intent="NIFTI_INTENT_NORMAL",
    )

    gii = nib.gifti.GiftiImage(darrays=[data_array])
    nib.save(
        gii,
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad.func.gii",
        ),
    )

    all_blurred = []
    # for i in ["midthickness-0mm", "0mm-1mm", "1mm-2mm", "2mm-3mm"]:
    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "-set-structure",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad.func.gii",
            ),
            "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
        ]
    )
    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "-metric-resample",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad.func.gii",
            ),
            # "E:\data\derivatives\micapipe\sub-PX103\ses-01\surf\sub-PX103_ses-01_Normal.func.gii",
            os.path.join(
                surf_dir,
                f"{bids_id}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii",
            ),
            f"src/data/fsLR-{resol}.{hemi}.sphere.reg.surf.gii",
            "BARYCENTRIC",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad-output.func.gii",
            ),
        ]
    )
    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "-set-structure",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad-output.func.gii",
            ),
            "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
        ]
    )
    return os.path.join(
        tmp_dir,
        f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad-output.func.gii",
    )


if __name__ == "__main__":
    sub = "sub-PX103"
    ses = "ses-01"
    surface = "fsnative"
    micapipe = "micapipe"
    hemi = "L"
    input_dir = f"E:/data/derivatives/{micapipe}/{sub}/{ses}/maps/"
    surf_dir = f"E:/data/derivatives/{micapipe}/{sub}/{ses}/surf/"
    output_dir = "."
    bids_id = f"{sub}_{ses}"
    compute_blurring(
        input_dir,
        surf_dir,
        bids_id,
        hemi,
        f"{output_dir}/{bids_id}_hemi-{hemi}_blurring.func.gii",
    )
