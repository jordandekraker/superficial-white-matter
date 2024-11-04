scans = [
    [
        "C:/Users/Ian/Documents/GitHub/Blurring/data/sub-PX001_ses-02_space-nativepro_T1w_brain.nii.gz",
        "E:/data/derivatives/micapipe/sub-PX001/ses-02/surf/sub-PX001_ses-02_hemi-R_space-nativepro_surf-fsLR-5k_label-",
    ],
    [
        "C:/Users/Ian/Documents/GitHub/Blurring/data/sub-PX006_ses-04_space-nativepro_T1w_brain.nii.gz",
        "E:/data/derivatives/micapipe/sub-PX006/ses-04/surf/sub-PX006_ses-04_hemi-L_space-nativepro_surf-fsLR-5k_label-",
    ],
    [
        "C:/Users/Ian/Documents/GitHub/Blurring/data/sub-PX029_ses-02_space-nativepro_T1w_brain.nii.gz",
        "E:/data/derivatives/micapipe/sub-PX029/ses-02/surf/sub-PX029_ses-02_hemi-L_space-nativepro_surf-fsLR-5k_label-",
    ],
    [
        "C:/Users/Ian/Documents/GitHub/Blurring/data/sub-PX049_ses-03_space-nativepro_T1w_brain.nii.gz",
        "E:/data/derivatives/micapipe/sub-PX049/ses-03/surf/sub-PX049_ses-03_hemi-L_space-nativepro_surf-fsLR-5k_label-",
    ],
    [
        "C:/Users/Ian/Documents/GitHub/Blurring/data/sub-PX051_ses-02_space-nativepro_T1w_brain.nii.gz",
        "E:/data/derivatives/micapipe/sub-PX051/ses-02/surf/sub-PX051_ses-02_hemi-L_space-nativepro_surf-fsLR-5k_label-",
    ],
    [
        "C:/Users/Ian/Documents/GitHub/Blurring/data/sub-PX069_ses-02_space-nativepro_T1w_brain.nii.gz",
        "E:/data/derivatives/micapipe/sub-PX069/ses-02/surf/sub-PX069_ses-02_hemi-L_space-nativepro_surf-fsLR-5k_label-",
    ],
]
import numpy as np
import os
import scipy
import subprocess
import nibabel as nib
import pandas as pd

wb_path = "C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64"
os.makedirs("output", exist_ok=True)
temppath = "output"
for scan in scans:
    micapath = scan[1]
    scan = scan[0]
    command = [
        "resseg-mni",
        scan,
        "-t",
        os.path.join(temppath, scan.split("/")[-1].split(".")[0] + "_resseg.tfm"),
    ]
    subprocess.run(command)
    command2 = [
        "resseg",
        scan,
        "-a",
        "20",
        "-t",
        os.path.join(temppath, scan.split("/")[-1].split(".")[0] + "_resseg.tfm"),
        "-o",
        os.path.join(temppath, scan.split("/")[-1].split(".")[0] + "_resseg.nii.gz"),
    ]
    subprocess.run(command2)

    command3 = [
        os.path.join(wb_path, "wb_command"),
        "-volume-to-surface-mapping",
        os.path.join(temppath, scan.split("/")[-1].split(".")[0] + "_resseg.nii.gz"),
        micapath + "midthickness.surf.gii",
        os.path.join(
            temppath, scan.split("/")[-1].split(".")[0] + "outputsurface.func.gii"
        ),
        "-trilinear",
    ]
    subprocess.run(command3)

    command4 = [
        os.path.join(wb_path, "wb_command"),
        "-metric-smoothing",
        "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/fsLR-5k.R.surf.gii",
        os.path.join(
            temppath, scan.split("/")[-1].split(".")[0] + "outputsurface.func.gii"
        ),
        "10",
        os.path.join(
            temppath, scan.split("/")[-1].split(".")[0] + "outputsurface.func.gii"
        ),
        "-fwhm",
        "-method",
        "GEO_GAUSS_EQUAL",
        "-fix-zeros",
    ]
    subprocess.run(command4)
    surface = (
        nib.load(
            os.path.join(
                temppath, scan.split("/")[-1].split(".")[0] + "outputsurface.func.gii"
            )
        )
        .darrays[0]
        .data
    )
    threshold = 0.01
    surface = np.where(surface < threshold, 0, 1)
    df = pd.DataFrame(surface)
    df.to_csv(
        os.path.join(
            temppath, scan.split("/")[-1].split(".")[0] + "thresholdedSurface.csv"
        ),
        index=False,
    )
    data_array = nib.gifti.gifti.GiftiDataArray(
        data=surface,
        intent="NIFTI_INTENT_NORMAL",
    )

    gii = nib.gifti.GiftiImage(darrays=[data_array])
    nib.save(
        gii,
        os.path.join(
            temppath, scan.split("/")[-1].split(".")[0] + "thresholdedSurface.func.gii"
        ),
    )
    print("e")
