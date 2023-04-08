from helper_json import loadJson
from helper_json import printJson


"""Select vertices based on their distances from the center of the colony (0,0).
If this distance is less than some predefined radius, add the vertex to the 
sub-colony being extracted. When a vertex is added to the subcolony, remove
any neighbors of its that are not also in the subcolony.

This is necessary in order to obtain a manageable colony size for gebc analysis."""


dataPath = "../../data/cellGrowth/simulation3/"
vertexCoord_fileName = dataPath + "vertexIdToCoords.json"
vertexIdToCoords = loadJson(vertexCoord_fileName)

voronoiSubNetwork_coords = {}
numVertices = 0
maxVertices = 5000

radius = 24
radiusSquare = radius**2

for id, coords in vertexIdToCoords.items():
    # if numVertices == maxVertices:
        # break
    x = float(coords[0])
    y = float(coords[1])
    if x**2 + y**2 < radiusSquare:
        voronoiSubNetwork_coords[id] = coords

print(f"Number of vertices in the extracted subnetwork: {len(voronoiSubNetwork_coords)}")

voronoiSubNetwork_coords_fileName = dataPath + "vertexIdToCoords_withinRadius"\
                                + str(radius) + ".json"
printJson(voronoiSubNetwork_coords, voronoiSubNetwork_coords_fileName)

vertexAdjList_fileName = dataPath + "vertexAdjList.json"
vertexAdjList = loadJson(vertexAdjList_fileName)
voronoiSubNetwork_adjList = {}

for id, neighbors in vertexAdjList.items():
    if id in voronoiSubNetwork_coords:
        # purge neighbors that are not in the subnetwork
        voronoiSubNetwork_adjList[id] = \
            [x for x in neighbors if str(x) in voronoiSubNetwork_coords]

voronoiSubNetwork_adjList_fileName = dataPath + "vertexAdjList_withinRadius"\
                                    + str(radius) + ".json"
printJson(voronoiSubNetwork_adjList, voronoiSubNetwork_adjList_fileName)
