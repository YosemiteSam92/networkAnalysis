# Computation of Geodesic Edge Betweenness Centrality for various types of network

Geodesic Edge Betweenness Centrality (gebc) is a measure of network centrality meant to indicate how crucial each vertex is to the flow of information in the network. A brief introduction to the concept can be found in the MS Power Point presentation in this folder. A good basic reference is also "Networks. An Introduction," by M. E. J. Newman. 

## Contents
*runScripts/epoxy*: 

- *analyseBetweennessCentrality.py*: computation of gebc for five networks made by cross-linking in-silico a DGEBA-DDS system with 127000 atoms (see https://pubs.acs.org/action/showCitFormats?doi=10.1021/acs.jctc.1c00423&ref=pdf). Input files in *data/epoxy/largestMolecularGroups_monomersOnly*; output in *results/epoxy*. 

- *analyseBetweennessCentrality_test.py*: computation of gebc for a few test networks, some of them with high geometrical symmetry.

- *reduceConnectivityToMonomers*:

    - *reduceConnectivityToMonomers.py*: Given the usual bond adjacency list containing all atoms, reduce its granularity by grouping all atom vertices belonging to a single monomer (either DGEBA or DDS) into a single vertex. See description in the script itself. Input files in *data/epoxy/lammpsData_annealedEquilibrated*; output in *data/epoxy/largestMolecularGroups_monomersOnly*. These reduced networks constitute the input for *analyseBetweennessCentrality.py*. 

    - *reduceConnectivityToMonomers_test.py*: same idea as the previous script; used to test the all-atom-to-monomer algorithm on a simple network. See *pdf* file for an explanation. 

- *other*: 
    - *removePbcsAndHfromDataFiles.py*: an earlier attempt to simplify the all-atom epoxy networks before computing gebc, by removing all H atoms and bonds straddling pbcs. The former turned out to be insufficient to make the problem tractable, while the latter turned out to be unnecessary for the purposes of computing gebc. 

*runScripts/cellGrowth*:

- *reconstructVoronoiConnectivity.py*: given tessellation data from Voronoi simulations of cell tissues (data/cellGrowth), reconstruct the network of Voronoi vertices. The output consists of three json files, with self-explanatory names, located in gebc/data/cellGrowth. In particular, the Voronoi connectivity is stored as adjacency list in *vertexAdjList.json*.

- *voronoiConnectivity.log*: the log output of *reconstructVoronoiConnectivity.py*, for debugging purposes.

- *analyseGebc.py*: calculate gebc values for the Voronoi vertices of the cell growth simulation. Input is *vertexAdjList.json*; output is *results/cellGrowth/voronoiVertices_gebc.dat*. 


*data*:

    input data for gebc calculations (see README in subfolders).

*results*:

    raw results of gebc calculations, in json format, associating a gebc score to each vertex id in the network.

- *test*: gebc analysis of various small networks, some of them with high symmetry, used for debugging purposes. 

*analysis*:

    analysis of raw gebc results, in the form of histograms and scatter plots.


*dot*:

    graph representations in dot format, for illustration/debugging purposes only.

*gebc*:

    the class that actually carries out the gebc calculation.
