import nibabel as nib
import numpy as np
from scipy.stats import mode
import subprocess
import os
from sWM import laplace_solver, surface_generator
import ants
from ants.ops import resample_image_to_target
import scipy
import pandas as pd


def fixmatrix(path, inputmap, outputmap, basemap, BIDS_ID, temppath, wb_path, mat_path):
    # Load the .mat file
    mat = scipy.io.loadmat(
        os.path.join(
            path,
            "xfm",
            f"{BIDS_ID}_{mat_path}.mat",
        )
    )

    # Extract variables from the .mat file
    affine_transform = mat["AffineTransform_double_3_3"].flatten()
    fixed = mat["fixed"].flatten()

    temp = np.identity(4)
    for i in range(3):
        temp[i, 3] = affine_transform[9 + i] + fixed[i]
        for j in range(3):
            temp[i, j] = affine_transform[i * 3 + j]
            temp[i, 3] -= temp[i, j] * fixed[j]

    flips = np.identity(4)
    flips[0, 0] = -1
    flips[1, 1] = -1

    m_matrix = np.linalg.inv(flips @ temp @ flips)

    print(m_matrix)
    with open(
        os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"),
        "w",
    ) as f:
        for row in m_matrix:
            f.write(" ".join(map(str, row)) + "\n")

    command4 = [
        os.path.join(wb_path, "wb_command"),
        "-convert-affine",
        "-from-world",
        os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"),
        "-to-world",
        "-inverse",
        os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"),
    ]

    subprocess.run(command4)

    command3 = [
        os.path.join(wb_path, "wb_command"),
        "-volume-resample",
        inputmap,
        basemap,
        "ENCLOSING_VOXEL",
        outputmap,
        "-affine",
        os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"),
    ]

    subprocess.run(command3)



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


def compute_blurring(
    input_dir,
    surf_dir,
    bids_id,
    hemi,
    feat,
    workbench_path,
    resol,
    fwhm,
    tmp_dir,
    fs_path,
    current_file_directory,
):

    base_path = input_dir
    for _ in range(4):  # Adjust the range to navigate up the desired number of levels
        base_path, _ = os.path.split(base_path)

    micapipe_path = os.path.split(input_dir)[0]

    freesurfer_path = os.path.join(
        base_path, "freesurfer", bids_id, "mri", "aparc+aseg.nii.gz"
    )

    temp_parc_path = os.path.join(
        tmp_dir, f"{bids_id}_{hemi}_surf-fsnative_label-temp.nii.gz"
    )
    print(temp_parc_path)
    output_path = os.path.join(tmp_dir, f"{bids_id}-laplace.nii.gz")
    # img = ants.image_read(freesurfer_path)
    # imgfixed = ants.image_read(f"{input_dir}/{bids_id}_space-nativepro_map-T1map.nii.gz")
    # resample_image_to_target(imgfixed, img, interp_type='multiLabel', verbose=True).to_filename(temp_parc_path)
    if not os.path.exists(freesurfer_path):
        subprocess.run(
            [
                os.path.join(fs_path, "mri_convert"),
                os.path.join(
                    micapipe_path,
                    f"{surf_dir}/mri/aparc+aseg.mgz",
                ),
                os.path.join(
                    tmp_dir, f"{bids_id}_{hemi}_surf-fsnative_label-temp-fixed.nii.gz"
                ),
            ]
        )
        freesurfer_path = os.path.join(
            tmp_dir, f"{bids_id}_{hemi}_surf-fsnative_label-temp-fixed.nii.gz"
        )
    if not os.path.exists(temp_parc_path):
        fixmatrix(
            path=input_dir,
            BIDS_ID=bids_id,
            temppath=tmp_dir,
            wb_path=workbench_path,
            inputmap=freesurfer_path,
            outputmap=temp_parc_path,
            basemap=f"{input_dir}/maps/{bids_id}_space-nativepro_map-T1map.nii.gz",
            mat_path="from-fsnative_to_nativepro_T1w_0GenericAffine",
        )

        if not os.path.exists(os.path.join(tmp_dir, "swm")):
            os.mkdir(os.path.join(tmp_dir, "swm"))
        laplace_solver.solve_laplace(temp_parc_path, output_path)
        surface_generator.shift_surface(
            f"{input_dir}/surf/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii",
            output_path,
            f"{tmp_dir}//swm//{bids_id}_{hemi}_sfwm-",
            [0.0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5],
        )

    wmBoundaryDataArr = load_gifti_data(
        f"{input_dir}/maps/{bids_id}_hemi-{hemi}_surf-fsnative_label-white_{feat}.func.gii"
    )
    wmBoundarySurfaceArr = load_gifti_data(
        f"{input_dir}/surf/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii"
    )

    modeofboundary = mode(wmBoundaryDataArr, keepdims=True)

    midthicknessDataArr = load_gifti_data(
        f"{input_dir}/maps/{bids_id}_hemi-{hemi}_surf-fsnative_label-midthickness_{feat}.func.gii"
    )
    midthicknessSurfaceArr = load_gifti_data(
        f"{input_dir}/surf/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-midthickness.surf.gii"
    )

    surfarr = [
        [midthicknessDataArr, midthicknessSurfaceArr],
        [wmBoundaryDataArr, wmBoundarySurfaceArr],
    ]
    if feat.lower() != "adc" or feat.lower() != "fa":
        volumemap = f"{input_dir}/maps/{bids_id}_space-nativepro_map-{feat}.nii.gz"
    else:
        volumemap = f"{input_dir}/maps/{bids_id}_space-nativepro_model-DTI_map-{feat.upper()}.nii.gz"
    for surf in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]:
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-volume-to-surface-mapping",
                volumemap,
                f"{tmp_dir}//swm//{bids_id}_{hemi}_sfwm-{surf}mm.surf.gii",
                f"{tmp_dir}//swm//{bids_id}_{hemi}_{feat}_{resol}_{fwhm}sfwm-{surf}mm-metric.func.gii",
                "-trilinear",
            ]
        )
        surfarr.append(
            [
                load_gifti_data(
                    f"{tmp_dir}//swm//{bids_id}_{hemi}_{feat}_{resol}_{fwhm}sfwm-{surf}mm-metric.func.gii"
                ),
                load_gifti_data(
                    f"{tmp_dir}//swm//{bids_id}_{hemi}_sfwm-{surf}mm.surf.gii"
                ),
            ]
        )

    distances = np.zeros(shape=(len(midthicknessDataArr), len(surfarr) - 1))
    dataArr = np.zeros(shape=(len(midthicknessDataArr), len(surfarr)))
    dataArr_nonmode = np.zeros(
        shape=(len(midthicknessDataArr), len(surfarr)), dtype=np.float32
    )
    for e, ds in enumerate(surfarr):
        data, surf = ds
        dataArr[:, e] = np.divide(data, modeofboundary.mode[0])
        dataArr_nonmode[:, e] = data
        if e == len(surfarr) - 1:
            break
        nextdata, nextsurt = surfarr[e + 1]
        print(e)
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

    data_non_grad = nib.gifti.gifti.GiftiDataArray(
        data=dataArr_nonmode,
        intent="NIFTI_INTENT_NORMAL",
    )

    gii_non_grad = nib.gifti.GiftiImage(darrays=[data_non_grad])
    nib.save(
        gii_non_grad,
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_NONgrad.func.gii",
        ),
    )

    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "-metric-resample",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_NONgrad.func.gii",
            ),
            os.path.join(
                input_dir,
                "surf",
                f"{bids_id}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii",
            ),
            os.path.join(
                current_file_directory,
                f"src/data/fsLR-{resol}.{hemi}.sphere.reg.surf.gii",
            ),
            "BARYCENTRIC",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_NONgrad-output.func.gii",
            ),
        ]
    )

    data_dist = nib.gifti.gifti.GiftiDataArray(
        data=distances.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL",
    )
    gii_dist = nib.gifti.GiftiImage(darrays=[data_dist])
    nib.save(
        gii_dist,
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_dist.func.gii",
        ),
    )

    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "-metric-resample",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_dist.func.gii",
            ),
            os.path.join(
                input_dir,
                "surf",
                f"{bids_id}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii",
            ),
            os.path.join(
                current_file_directory,
                f"src/data/fsLR-{resol}.{hemi}.sphere.reg.surf.gii",
            ),
            "BARYCENTRIC",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_dist-output.func.gii",
            ),
        ]
    )

    data_fslr = nib.load(
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_NONgrad-output.func.gii",
        )
    ).darrays
    data_fslr = [x.data for x in data_fslr]
    print(data_fslr)

    df = pd.DataFrame(data_fslr)
    df.to_csv(
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-intensities.csv",
        ),
        index=False,
    )

    data_dist = nib.load(
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_dist-output.func.gii",
        )
    ).darrays
    data_dist = [x.data for x in data_dist]

    distancesdf = pd.DataFrame(data_dist)
    distancesdf.to_csv(
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-distances.csv",
        ),
        index=False,
    )

    return [
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_NONgrad-output.func.gii",
        ),
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-intensities.csv",
        ),
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-distances.csv",
        ),
    ]


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
