from helper_json import loadJson
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

"""Extract and plot the distribution of the number of 
neighbors of each vertex"""

dataPath = "../../../data/cellGrowth/"
radius = "30"
adjList_fileName = dataPath + "vertexAdjList_withinRadius"\
    + radius + ".json"
adjList = loadJson(adjList_fileName)

distOfNeighbors = [
    len(lisOfNeighs) for lisOfNeighs in adjList.values()
]
distOfNeighbors = np.array(distOfNeighbors, dtype=int)

fig, ax = plt.subplots()
nBins = 80
ax.set_xlabel('number of neighbors of each Voronoi vertex')
ax.set_ylabel('frequency')
histValues, bins, patches = ax.hist(
    distOfNeighbors, bins=nBins, density=False)

fileName = "numOfNeighbors_hist_withinRadius" + radius
plt.savefig(fileName + ".pdf", format="pdf", dpi=300)

