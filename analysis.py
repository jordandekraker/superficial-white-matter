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


hemis = ["L", "R"]
datadir = "E:\data\derivatives\zblur"
for fold in os.listdir(datadir):
    delete_empty_folder(os.path.join(datadir, fold))

totalsessions = []
for fold in os.listdir(datadir):
    sessions = os.listdir(os.path.join(datadir, fold))
    output_sessions = []
    for session in sessions:
        start = session.find(f"{fold}_") + len(f"{fold}_")
        end = session.find("_L") if "_L" in session else session.find("_R")
        session_number = session[start:end]
        print(f"Session: {session}, Extracted Substring: {session_number}")
        output_sessions.append(session_number)
    sessions = list(set(output_sessions))
    for ses in sessions:
        totalsessions.append(ses)


totalsubs = len(totalsessions)
L_intensities_array = np.zeros((4842, 12, totalsubs))
L_distances_array = np.zeros((4842, 11, totalsubs))

R_intensities_array = np.zeros((4842, 12, totalsubs))
R_distances_array = np.zeros((4842, 11, totalsubs))

e = 0
for fold in os.listdir(datadir):
    sessions = os.listdir(os.path.join(datadir, fold))
    output_sessions = []
    for session in sessions:
        start = session.find(f"{fold}_") + len(f"{fold}_")
        end = session.find("_L") if "_L" in session else session.find("_R")
        session_number = session[start:end]
        print(f"Session: {session}, Extracted Substring: {session_number}")
        output_sessions.append(session_number)
    sessions = list(set(output_sessions))
    for session in sessions:
        for hemi in hemis:
            distances_path = os.path.join(
                datadir,
                fold,
                f"sub-{fold}_{session}_{hemi}_T1map_blur_distances.csv",
            )
            intensities_path = os.path.join(
                datadir,
                fold,
                f"sub-{fold}_{session}_{hemi}_T1map_blur_intensities.csv",
            )
            distances = np.genfromtxt(
                distances_path, delimiter=",", skip_header=1
            ).transpose()
            intensities = np.genfromtxt(
                intensities_path, delimiter=",", skip_header=1
            ).transpose()
            if hemi == "L":
                L_intensities_array[:, :, int(e)] = intensities
                L_distances_array[:, :, int(e)] = distances
            else:
                R_intensities_array[:, :, int(e)] = intensities
                R_distances_array[:, :, int(e)] = distances
        e += 1

LIavgacrosstrial = np.mean(L_intensities_array, axis=0).transpose()
LIstdacrosstrial = np.std(L_intensities_array, axis=0).transpose()

LDavgacrosstrial = np.mean(L_distances_array, axis=0).transpose()
LDstdacrosstrial = np.std(L_distances_array, axis=0).transpose()


LDavgacrosstrial_reshaped = np.zeros((len(LDavgacrosstrial), 12))
for en, x in enumerate(LDavgacrosstrial):
    for e in range(len(x)):
        if e == 0:
            LDavgacrosstrial_reshaped[en, e] = -x[e]
            LDavgacrosstrial_reshaped[en, e + 1] = 0
        else:
            LDavgacrosstrial_reshaped[en, e + 1] = (
                LDavgacrosstrial_reshaped[en, e] + x[e]
            )
    print(x)


threshold = 5.5
mask = (
    LDavgacrosstrial_reshaped[:, 11] <= threshold
)  # Assuming you want to check the first trial


LIavgacrosstrial = LIavgacrosstrial[mask]
LIstdacrosstrial = LIstdacrosstrial[mask]

LDavgacrosstrial = LDavgacrosstrial[mask]
LDstdacrosstrial = LDstdacrosstrial[mask]
LDavgacrosstrial_reshaped = LDavgacrosstrial_reshaped[mask]

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

for row, std_dev, dists in zip(
    LIavgacrosstrial, LIstdacrosstrial, LDavgacrosstrial_reshaped
):
    plt.errorbar(dists, row, yerr=std_dev, marker="o", capsize=5)

plt.xlabel("Distance")
plt.ylabel("Value")
plt.title("Line Graph for Each Row in LIavgacrosstrial with Standard Deviations")
# plt.legend([f"Row {i+1}" for i in range(LIavgacrosstrial.shape[0])], loc="upper right")
plt.grid(True)
plt.show()


LIavgacrosstrial = np.mean(L_intensities_array, axis=2)
LIstdacrosstrial = np.std(L_intensities_array, axis=2)

LDavgacrosstrial = np.mean(L_distances_array, axis=2)
LDstdacrosstrial = np.std(L_distances_array, axis=2)

LDavgacrosstrial_reshaped = np.zeros((len(LDavgacrosstrial), 12))
for en, x in enumerate(LDavgacrosstrial):
    for e in range(len(x)):
        if e == 0:
            LDavgacrosstrial_reshaped[en, e] = -x[e]
            LDavgacrosstrial_reshaped[en, e + 1] = 0
        else:
            LDavgacrosstrial_reshaped[en, e + 1] = (
                LDavgacrosstrial_reshaped[en, e] + x[e]
            )
    print(x)
threshold = 5.5
mask = (
    LDavgacrosstrial_reshaped[:, 11] <= threshold
)  # Assuming you want to check the first trial


LIavgacrosstrial = LIavgacrosstrial[mask]
LIstdacrosstrial = LIstdacrosstrial[mask]

LDavgacrosstrial = LDavgacrosstrial[mask]
LDstdacrosstrial = LDstdacrosstrial[mask]
LDavgacrosstrial_reshaped = LDavgacrosstrial_reshaped[mask]
# test = LDavgacrosstrial_reshaped[-500:-1]
# test2 = LDavgacrosstrial_reshaped[-1000:-501]
# test3 = LIavgacrosstrial[-500:-1]
# test4 = LIavgacrosstrial[-1000:-501]


import matplotlib.pyplot as plt


# Plotting
plt.figure(figsize=(10, 6))

# Set x-axis limits
plt.xlim(-3, 6)

for row, std_dev, dists in zip(
    LIavgacrosstrial, LIstdacrosstrial, LDavgacrosstrial_reshaped
):
    plt.errorbar(dists, row, yerr=std_dev, marker="o", capsize=5)
# LDavgacrosstrial_reshaped = LDavgacrosstrial_reshaped[mask]
plt.xlabel("Distance")
plt.ylabel("Value")
plt.title("Line Graph for Each Row in LIavgacrosstrial with Standard Deviations")
# plt.legend([f"Row {i+1}" for i in range(LIavgacrosstrial.shape[0])], loc="upper right")
plt.grid(True)
plt.show()
print("e")
