# Analysis of tessellation data from Voronoi simulations of cell tissues 

*voronoi_cluster_simulations.odp* summarizes the simulation results for many cell clusters. The simulation numbers naming the subfolders here correspond to the simulation number in this file.

We want to compare the gebc behavior of colonies showing different distributions of neighbors for the cells in the bulk of the cluster. In particular, simulation 4 shows the typical spike at 6 neighbors, while simulation 3 also presents a noteworthy spike at 5. This same atypical peak is also found in other systems. Could this be related to information flow in the network? 

Gebc can assess this at the level of geodesic communication among Voronoi vertices.

## Folders

- pilotSimulation: 

    the first batch of Voronoi data analyzed in the summer of 2022 while developing the workflow. Files:

    *neighbors...txt*:
    
    cell_id cell_id_of_neighbor_1 cell_id_of_neighbor_2 ...

    *voronoi...txt*:
    
    cell_id number_of_neighbors x_coord._of_vertex_1 y_coord._of_vertex 1 x_coord._of_vertex_2 y_coord._of_vertex_2...

    Vertices that are adjacent in this list are connected; the last vertex on the right is identical to the first on the left.


    *nuclei...txt*:

    Information about cell nuclei. Unnecessary to compute the gebc of the vertices.

    *cellIdToVertexIds.json*, *vertexIdToCoords.json*:

    The output of reconstructVoronoiConnectivity.py in runScripts. The names are self-explanatory.

    *vertexIdToCoords_withinRadiusn.json, vertexAdjList_withinRadiusn.json*:

    Same as above, but for a subcolony centered at the origin, bounded by a certain radius *n*. Extracting a subnetwork was necessary because the whole colony was too large for gebc calculation. 

- simulation3 and simulation4: 
    same file structure and file formats as pilotSimulation. 



