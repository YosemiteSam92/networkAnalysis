from os import chdir, path
from helper_json import loadJson
import numpy as np
import matplotlib.pyplot as plt


# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

gebcs_fileNames = ["largestGroup_monomers_large" + str(i) for i in range(1, 6)]
gebcs = [loadJson("../../../results/epoxy/" + fileName + "_gebc.json")
         for fileName in gebcs_fileNames]

nBins = 80

# unified gebc data spanning all systems
gebc_allSys = []

# individual system histograms, all in one file

index = 0
fig, axes = plt.subplots(3, 2)

for axRow in axes:
    for ax in axRow:

        gebc = gebcs[index]
        values_list = list(gebc.values())
        gebc_allSys += values_list
        values = np.array(values_list)
        values /= max(values)

        ax.set_title(str(index + 1), y=.7)
        histValues, bins, _ = ax.hist(values, bins=nBins,
                                      density=True)

        fileName = "normalizedGebc_histDensity_large" + str((index + 1))
        np.savetxt(fileName + ".txt", histValues)

        index += 1
        if index > len(gebcs) - 1:
            break


for ax in axes.flat:
    # keep x labels only for lower graphs
    # and y labels only for left graphs
    ax.label_outer()

# suppress lower-right graph
axes[2, 1].set_axis_off()

fig.supxlabel("normalized gebc")
fig.supylabel("probability density")

fileName = "normalizedGebc_histDensity_largeAll"
plt.savefig(fileName + ".pdf", format="pdf", dpi=300)


# unified data histogram

gebc_allSys = np.array(gebc_allSys)
gebc_allSys /= max(gebc_allSys)

fig, ax = plt.subplots(1, 1)
fig.suptitle("All systems aggregated, " + str(nBins) + " bins", y=0.95)
ax.set_ylabel("probability density")
ax.set_xlabel("normalized GEBC")

hist, bin_edges, _ = plt.hist(gebc_allSys, bins=nBins, density=True)
fileName = "normalizedGebc_histDensity_large_AllSysAggregated"
np.savetxt(fileName + ".txt", hist)
plt.savefig(fileName + ".pdf", format="pdf", dpi=300)
