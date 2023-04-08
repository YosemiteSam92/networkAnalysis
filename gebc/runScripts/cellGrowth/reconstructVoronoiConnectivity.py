import sys
from os.path import dirname, abspath
import argparse
from helper import getStripSplitLine
from helper_json import printJson
from helper_dict import getFirstKeyFromValue
from helper_dict import changeTypeOfDictKeys

sys.path.append(abspath(dirname(__file__)))

parser = argparse.ArgumentParser(
        description="Reconstruct the network of Voronoi vertices of the cluster."
    )
parser.add_argument(
    "-d",
    "--dataPath",
    default="../../data/cellGrowth/",
    help="The directory for input and output files.",
)
parser.add_argument(
    "-n",
    "--neighbors",
    default="",
    help="File name of input data about neighboring cells, without extension."
         + "It must be located in dataPath.",
)
parser.add_argument(
    "-v",
    "--voronoi",
    default="",
    help="File name of input data about Voronoi vertices, without extension."
         + "It must be located in dataPath.",
)

args = parser.parse_args()
dataPath = args.dataPath
neighborsFileName = args.neighbors
voronoiFileName = args.voronoi


"""read table of neighbors
   different rows have different numbers of entries,
   so np.loadtxt will not work"""

# cellIdToNeighborsIds = {}

# with open(dataPath + neighborsFileName + ".txt", "r") as file:
#     while 1:
#         lineList = getStripSplitLine(file)
#         if lineList == [""]:
#             break
#         cellIdToNeighborsIds[lineList[0]] = lineList[1:]

"""read table of neighbors' coordinates"""

cellIdToNeighborsCoords = {}

with open(dataPath + voronoiFileName + ".txt", "r") as file:
    while 1:
        lineList = getStripSplitLine(file)
        if lineList == [""]:
            break
        cellIdToNeighborsCoords[lineList[0]] = lineList[1:]

"""associate each cellId with its list of vertexIds
   associate each vertexId with the vertex's coord.s
   find the neighbors of each vertex"""

vertexIdToCoords = {}
cellIdtoVertexIds = {}
vertexAdjList = {}
nextNewVertexId = 0
prevId = 0
newVertex = False

for cellId, neighCoords in cellIdToNeighborsCoords.items():
    # if cellId == "52003":
    #     break
    cellIdtoVertexIds[cellId] = []
    numVertices = int(neighCoords[0]) - 1  # last pair is identical to first
    print("cell id: ", cellId)
    # print("num. vertices: ", numVertices)
    count = 0
    while count < (numVertices*2):

        vertexCoords = (neighCoords[1+count], neighCoords[2+count])

        # is this a new vertex? check if coord.s not seen before
        # this search is the bottleneck
        if vertexCoords not in vertexIdToCoords.values():
            newVertex = True
            vertexId = nextNewVertexId
            vertexIdToCoords[vertexId] = vertexCoords
            vertexAdjList[vertexId] = []
            # print("this vertex is new:", vertexId)
        else:
            vertexId = getFirstKeyFromValue(vertexIdToCoords, vertexCoords)
            # print("this vertex is known:", vertexId)

        cellIdtoVertexIds[cellId].append(vertexId)

        # update connectivity of this vertex (whether new or not)
        # print("count: ", count)
        if count > 0:
            if prevId not in vertexAdjList[vertexId]:
                vertexAdjList[vertexId].append(prevId)
                vertexAdjList[prevId].append(vertexId)
        count += 2

        # print(f"neighbors of {vertexId}: {vertexAdjList[vertexId]}")
        # print(f"neighbors of {prevId}: {vertexAdjList[prevId]}")

        prevId = vertexId

        if newVertex:
            nextNewVertexId += 1
            newVertex = False

        # print("cell vertices:", cellIdtoVertexIds[cellId])
        # print(f"coord.s of {vertexId}:", vertexIdToCoords[vertexId])

    # last vertex is linked to the first one
    firstVertexOfCell = cellIdtoVertexIds[cellId][0]
    if firstVertexOfCell not in vertexAdjList[vertexId]:
        vertexAdjList[vertexId].append(firstVertexOfCell)
    # print(f"neighbors of {vertexId}: {vertexAdjList[vertexId]}")
    # print(f"all coords of vertices: {printList(list(vertexIdToCoords.values()))}")
    # break

"""print dictionaries to file"""

printJson(changeTypeOfDictKeys(cellIdtoVertexIds, str, int),
          dataPath + "cellIdToVertexIds.json")
printJson(vertexIdToCoords, dataPath + "vertexIdToCoords.json")
printJson(changeTypeOfDictKeys(vertexAdjList, str, int),
          dataPath + "vertexAdjList.json")
