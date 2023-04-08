from helper_json import loadJson
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

resultsPath = "../../../results/cellGrowth/"
fileName = "voronoiVertices_withinRadius30_gebc"
vertices_gebc = loadJson(resultsPath + fileName + ".json")

gebc = np.array(list(vertices_gebc.values()))
gebc /= max(gebc)

fig, ax = plt.subplots()
nBins = 80
ax.set_xlabel('normalized GEBC')
ax.set_ylabel('probability density')
histValues, bins, patches = ax.hist(gebc, bins=nBins, density=True)

# color code the histogram
fracs = histValues/histValues.max()
norm = colors.Normalize(fracs.min(), fracs.max())
for frac, patch in zip(fracs, patches):
    color = plt.cm.viridis_r(norm(frac))
    patch.set_facecolor(color)

fileName = "normalizedGebc_histDensity_withinRadius30"
plt.savefig(fileName + ".pdf", format="pdf", dpi=300)
