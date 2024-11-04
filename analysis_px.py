import pandas as pd
import os
import numpy as np


def delete_empty_folder(folder_path):
    # Check if the folder is empty
    if not os.listdir(folder_path):
        # Delete the folder if it is empty
        os.rmdir(folder_path)
        print(f"Deleted empty folder: {folder_path}")
    else:
        print(f"Folder is not empty: {folder_path}")


L_intensities_array = np.zeros((4842, 12))
L_distances_array = np.zeros((4842, 11))

R_intensities_array = np.zeros((4842, 12))
R_distances_array = np.zeros((4842, 11))

distances = np.genfromtxt(
    "E:/data/derivatives/zblur_patients/PX001/sub-PX001_ses-01_L_T1map_blur_distances.csv",
    delimiter=",",
    skip_header=1,
).transpose()
intensities = np.genfromtxt(
    "E:/data/derivatives/zblur_patients/PX001/sub-PX001_ses-01_L_T1map_blur_intensities.csv",
    delimiter=",",
    skip_header=1,
).transpose()

mask = np.genfromtxt(
    "C:/Users/Ian/Documents/GitHub/Blurring/output/sub-PX001_ses-02_space-nativepro_T1w_brainthresholdedSurface.csv",
    delimiter=",",
    skip_header=1,
).transpose()

# Convert mask to boolean array
mask = mask.astype(bool)

R_intensities_array[:, :] = intensities
R_distances_array[:, :] = distances


RIavgacrosstrial = R_intensities_array
RDavgacrosstrial = R_distances_array


RDavgacrosstrial_reshaped = np.zeros((len(RDavgacrosstrial), 12))
for en, x in enumerate(RDavgacrosstrial):
    for e in range(len(x)):
        if e == 0:
            RDavgacrosstrial_reshaped[en, e] = -x[e]
            RDavgacrosstrial_reshaped[en, e + 1] = 0
        else:
            RDavgacrosstrial_reshaped[en, e + 1] = (
                RDavgacrosstrial_reshaped[en, e] + x[e]
            )
    print(x)


# threshold = 5.5
# mask = (
#     LDavgacrosstrial_reshaped[:, 11] <= threshold
# )  # Assuming you want to check the first trial


RIavgacrosstrial = RIavgacrosstrial[mask]


RDavgacrosstrial = RDavgacrosstrial[mask]

RDavgacrosstrial_reshaped = RDavgacrosstrial_reshaped[mask]

# test = LDavgacrosstrial_reshaped[-500:-1]
# test2 = LDavgacrosstrial_reshaped[-1000:-501]
# test3 = LIavgacrosstrial[-500:-1]
# test4 = LIavgacrosstrial[-1000:-501]


import matplotlib.pyplot as plt

plt.close()

# Plotting
plt.figure(figsize=(10, 6))
plt.clf()
# Set x-axis limits
plt.xlim(-3, 6)

for row, dists in zip(RIavgacrosstrial, RDavgacrosstrial_reshaped):
    plt.errorbar(dists, row, marker="o", capsize=5)

plt.xlabel("Distance")
plt.ylabel("Value")
plt.title("Line Graph for Each Row in LIavgacrosstrial with Standard Deviations")
# plt.legend([f"Row {i+1}" for i in range(LIavgacrosstrial.shape[0])], loc="upper right")
plt.grid(True)
plt.show()


# LIavgacrosstrial = np.mean(L_intensities_array, axis=2)
# LIstdacrosstrial = np.std(L_intensities_array, axis=2)

# LDavgacrosstrial = np.mean(L_distances_array, axis=2)
# LDstdacrosstrial = np.std(L_distances_array, axis=2)

# LDavgacrosstrial_reshaped = np.zeros((len(LDavgacrosstrial), 12))
# for en, x in enumerate(LDavgacrosstrial):
#     for e in range(len(x)):
#         if e == 0:
#             LDavgacrosstrial_reshaped[en, e] = -x[e]
#             LDavgacrosstrial_reshaped[en, e + 1] = 0
#         else:
#             LDavgacrosstrial_reshaped[en, e + 1] = (
#                 LDavgacrosstrial_reshaped[en, e] + x[e]
#             )
#     print(x)
# threshold = 5.5
# mask = (
#     LDavgacrosstrial_reshaped[:, 11] <= threshold
# )  # Assuming you want to check the first trial


# LIavgacrosstrial = LIavgacrosstrial[mask]
# LIstdacrosstrial = LIstdacrosstrial[mask]

# LDavgacrosstrial = LDavgacrosstrial[mask]
# LDstdacrosstrial = LDstdacrosstrial[mask]
# LDavgacrosstrial_reshaped = LDavgacrosstrial_reshaped[mask]
# # test = LDavgacrosstrial_reshaped[-500:-1]
# # test2 = LDavgacrosstrial_reshaped[-1000:-501]
# # test3 = LIavgacrosstrial[-500:-1]
# # test4 = LIavgacrosstrial[-1000:-501]


# import matplotlib.pyplot as plt


# # Plotting
# plt.figure(figsize=(10, 6))

# # Set x-axis limits
# plt.xlim(-3, 6)

# for row, std_dev, dists in zip(
#     LIavgacrosstrial, LIstdacrosstrial, LDavgacrosstrial_reshaped
# ):
#     plt.errorbar(dists, row, yerr=std_dev, marker="o", capsize=5)
# # LDavgacrosstrial_reshaped = LDavgacrosstrial_reshaped[mask]
# plt.xlabel("Distance")
# plt.ylabel("Value")
# plt.title("Line Graph for Each Row in LIavgacrosstrial with Standard Deviations")
# # plt.legend([f"Row {i+1}" for i in range(LIavgacrosstrial.shape[0])], loc="upper right")
# plt.grid(True)
# plt.show()
# print("e")
