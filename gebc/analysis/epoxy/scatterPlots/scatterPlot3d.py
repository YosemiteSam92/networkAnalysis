from lammpsData import LammpsData
from helper_json import loadJson
import numpy as np
import matplotlib.pyplot as plt

gebcs_fileNames = ["largestGroup_monomers_large" + str(i) for i in range(1, 6)]

dataFiles = [
    "data.equilibriate_298K_large1_127ns",
    "data.equilibriate_298K_large2_152ns",
    "data.equilibriate_298K_large3_152ns",
    "data.equilibriate_298K_large4_152ns",
    "data.equilibriate_298K_large5_107ns"
]

# individual scatter plots for each system

# for sysNum in range(len(dataFiles)):

#     data = LammpsData("../../../data/epoxy/lammpsData_annealedEquilibrated/"
#                       + dataFiles[sysNum])
#     gebc = loadJson("../../../results/epoxy/" + gebcs_fileNames[sysNum]
#                     + "_gebc.json")

#     # get coordinates of monomers (S for DDS and central C for DGEBA)
#     representativeAtomIdsOfMonomers = list(gebc.keys())  # strings
#     representativeAtomIndicesOfMonomers = np.array(
#         representativeAtomIdsOfMonomers).astype(int) - 1
#     monomersX = data.x[representativeAtomIndicesOfMonomers]
#     # slicing respects the order of the indices in the array of indices

#     # normalize gebcs values before passing them as color map
#     gebcValues = np.array(list(gebc.values()))
#     gebcValues /= gebcValues.max()

#     fig = plt.figure(figsize=(12, 9.6))
#     ax = fig.add_subplot(projection='3d')
#     plot = ax.scatter(monomersX[:, 0], monomersX[:, 1], monomersX[:, 2],
#                       marker="o",
#                       c=gebcValues,
#                       cmap="viridis_r")
#     # _r reverses the color map, in order to exalt high gebc values better
#     fig.colorbar(plot, shrink=0.7)
#     ax.set_xlabel('X / Å')
#     ax.set_ylabel('Y / Å')
#     ax.set_zlabel('Z / Å')
#     plt.savefig(gebcs_fileNames[sysNum] + "_gebc.pdf", format="pdf", dpi=300)

# All scatter plots in 1 figure

sysNum = 0

fig = plt.figure(figsize=(20, 16))

# get rid of white margins (not really working)
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0

labelsFont = 14

for sysNum in range(len(dataFiles)):

    ax = fig.add_subplot(3, 2, sysNum + 1, projection='3d')
    data = LammpsData("../../../data/epoxy/lammpsData_annealedEquilibrated/"
                      + dataFiles[sysNum])
    gebc = loadJson("../../../results/epoxy/"
                    + gebcs_fileNames[sysNum] + "_gebc.json")

    # get coordinates of monomers
    # (representative atoms are S for DDS and central C for DGEBA)
    representativeAtomIdsOfMonomers = list(gebc.keys())  # strings
    # cast to np array of int
    representativeAtomIndicesOfMonomers = np.array(
        representativeAtomIdsOfMonomers).astype(int) - 1
    # slice to retain only the coord.s of rep. atoms
    # (slicing respects the order of the indices in the array of indices)
    monomersX = data.x[representativeAtomIndicesOfMonomers]

    # normalize gebcs values before passing them as color map
    gebcValues = np.array(list(gebc.values()))
    gebcValues /= gebcValues.max()

    plot = ax.scatter(monomersX[:, 0], monomersX[:, 1], monomersX[:, 2],
                      marker="o",
                      c=gebcValues,
                      cmap="viridis_r")

    ax.set_title(str(sysNum + 1), y=1, fontsize=20)
    ax.set_xlabel('X / Å', fontsize=labelsFont)
    ax.set_ylabel('Y / Å', fontsize=labelsFont)
    ax.set_zlabel('Z / Å', fontsize=labelsFont)
    ax.tick_params(axis='x', labelsize=labelsFont-2)
    ax.tick_params(axis='y', labelsize=labelsFont-2)
    ax.tick_params(axis='z', labelsize=labelsFont-2)

    if sysNum > len(dataFiles) - 1:
        break

cbar = fig.colorbar(plot, shrink=0.7, pad=0.1)
cbar.set_label("gebc", fontsize=labelsFont+5)
cbar.ax.tick_params(labelsize=labelsFont)
# _r reverses the color map, in order to exalt high gebc values better

fig.tight_layout()  # automatically remove unused white space
fileName = "largestGroup_monomers_gebcColor_scatter3d_largeAll"
plt.savefig(fileName + ".pdf", format="pdf", dpi=300)
