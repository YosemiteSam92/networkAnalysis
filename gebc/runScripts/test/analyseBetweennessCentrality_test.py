from geodesicEdgeBetweenneessCentrality\
    import GeodesicEdgeBetweennessCentrality as gebc
from os import chdir, path

# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

outFolder = "../../results/test/"

# test graphs
# g = {}
# g[0] = [1, 2]
# g[1] = [0, 3, 4]
# g[2] = [0, 3, 4]
# g[3] = [2, 1, 5]
# g[4] = [2, 1, 5]
# g[5] = [3, 4]
# gebcAnalyzer = gebc(g, parallel=False)
# expected gebc (by hand):
# [0.33333333 1.83333333 1.83333333 1.83333333 1.83333333 0.33333333]

# g = {}
# g[0] = [1, 4, 6]
# g[1] = [0, 2]
# g[2] = [1, 3]
# g[3] = [2, 4]
# g[4] = [3, 5, 0]
# g[5] = [4, 6]
# g[6] = [5, 0]
# gebcAnalyzer = gebc(g, parallel=False)
# expected gebc (by hand):
# [5. 2. 1. 2. 5. 1. 1.]

# g = {}
# g[0] = [1, 2, 3]
# g[1] = [0, 2, 4]
# g[2] = [0, 1, 3, 4]
# g[3] = [0, 2, 4]
# g[4] = [1, 2, 3, 5]
# g[5] = [4, 7]
# g[6] = [7]
# g[7] = [5, 6, 8]
# g[8] = [7]
# gebcAnalyzer = gebc(g, parallel=False)

# In the following tests, check that geometrical symmetries
# are reflected in the set of gebc values

# two adjacent rhombi
# g = {}
# g[0] = [1, 3]
# g[1] = [0, 2]
# g[2] = [1, 3, 4, 6]
# g[3] = [0, 2]
# g[4] = [2, 5]
# g[5] = [4, 6]
# g[6] = [2, 5]
# gebcAnalyzer = gebc(g, parallel=False)
# gebcAnalyzer.print(outFolder + "twoRhombi_gebc.txt")

# hexagon
# g = {}
# g[0] = [1, 5]
# g[1] = [0, 2]
# g[2] = [1, 3]
# g[3] = [2, 4]
# g[4] = [3, 5]
# g[5] = [4, 0]
# gebcAnalyzer = gebc(g, parallel=False)
# gebcAnalyzer.print(outFolder + "hexagon_gebc.txt")

# two adjacent hexagons
# g = {}
# g[0] = [1, 5]
# g[1] = [0, 2]
# g[2] = [1, 3, 6]
# g[3] = [2, 4, 9]
# g[4] = [3, 5]
# g[5] = [4, 0]
# g[6] = [2, 7]
# g[7] = [6, 8]
# g[8] = [7, 9]
# g[9] = [8, 3]
# gebcAnalyzer = gebc(g, parallel=False)
# gebcAnalyzer.print(outFolder + "twoHexagons_gebc.txt")

# three adjacent hexagons
# g = {}
# g[0] = [1, 5]
# g[1] = [0, 2]
# g[2] = [1, 3]
# g[3] = [2, 4, 12]
# g[4] = [3, 5, 9]
# g[5] = [0, 4, 6]
# g[6] = [5, 7]
# g[7] = [6, 8]
# g[8] = [7, 9]
# g[9] = [4, 8, 10]
# g[10] = [9, 11]
# g[11] = [10, 12]
# g[12] = [3, 11]
# gebcAnalyzer = gebc(g, parallel=False)
# gebcAnalyzer.print(outFolder + "threeHexagons_gebc.txt")

# four adjacent hexagons
# g = {}
# g[0] = [1, 5]
# g[1] = [0, 2]
# g[2] = [1, 3]
# g[3] = [2, 4, 12]
# g[4] = [3, 5, 9]
# g[5] = [0, 4, 6]
# g[6] = [5, 7]
# g[7] = [6, 8]
# g[8] = [7, 9, 15]
# g[9] = [4, 8, 10]
# g[10] = [9, 11, 13]
# g[11] = [10, 12]
# g[12] = [3, 11]
# g[13] = [10, 14]
# g[14] = [13, 15]
# g[15] = [8, 14]
# gebcAnalyzer = gebc(g, parallel=False)
# gebcAnalyzer.print(outFolder + "fourHexagons_gebc.txt")

# three adjacent hexagons and one pentagon ("5-point" defect)
# g = {}
# g[0] = [1, 5]
# g[1] = [0, 2]
# g[2] = [1, 3]
# g[3] = [2, 4, 11]
# g[4] = [3, 5, 8]
# g[5] = [0, 4, 6]
# g[6] = [5, 7]
# g[7] = [6, 8, 14]
# g[8] = [4, 7, 9]
# g[9] = [8, 10, 12]
# g[10] = [9, 11]
# g[11] = [3, 10]
# g[12] = [9, 13]
# g[13] = [12, 14]
# g[14] = [7, 13]
# gebcAnalyzer = gebc(g, parallel=False)
# gebcAnalyzer.print(outFolder + "fourHexagons_5pointDefect_gebc.txt")

# three adjacent hexagons and one heptagon ("7-point" defect)
g = {}
g[0] = [1, 5]
g[1] = [0, 2]
g[2] = [1, 3]
g[3] = [2, 4, 13]
g[4] = [3, 5, 10]
g[5] = [0, 4, 6]
g[6] = [5, 7]
g[7] = [6, 8]
g[8] = [7, 9]
g[9] = [8, 10, 16]
g[10] = [4, 9, 11]
g[11] = [10, 12, 14]
g[12] = [11, 13]
g[13] = [3, 12]
g[14] = [11, 15]
g[15] = [14, 16]
g[16] = [9, 15]
gebcAnalyzer = gebc(g, parallel=False)
gebcAnalyzer.print(outFolder + "fourHexagons_7pointDefect_gebc.txt")
