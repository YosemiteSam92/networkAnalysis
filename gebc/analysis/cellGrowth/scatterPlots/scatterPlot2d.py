from helper_json import loadJson
import matplotlib.pyplot as plt

resultsPath = "../../../results/cellGrowth/simulation3/"
dataPath = "../../../data/cellGrowth/simulation3/"
plotPath = "simulation3/"

radius = "37"
fileName = "voronoiVertices_withinRadius" + radius + "_unnormalized_gebc"
vertices_gebc = loadJson(resultsPath + fileName + ".json")

fileName = "vertexIdToCoords_withinRadius" + radius
vertices_coords = loadJson(dataPath + fileName + ".json")

gebc = list(vertices_gebc.values())
coords = list(vertices_coords.values())

x = [float(coord[0]) for coord in coords]
y = [float(coord[1]) for coord in coords]

fig, ax = plt.subplots()
plot = ax.scatter(x, y, marker="o", c=gebc, cmap="viridis_r")
fig.colorbar(plot, shrink=0.7)
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.savefig(plotPath + "voronoiVertices_withinRadius" + radius
            + "_unnormalized_gebc.pdf", format="pdf", dpi=300)
