'''
gebc(vertex v) = sum_(i<j) n_(ij)^v / g_ij

where

    n_(ij)^v is defined to be 1 if there exists a geodesic
    path connecting vertices i and j passing through v and 0 otherwise

    g_ij is the total number of geodesic paths between i and j

so each path has a weight equal to the inverse of the
number of paths.

Useful resources on BFS and reconstructing all (not just one)
geodesic paths from it:
https://www.youtube.com/watch?v=09_LlHjoEiY (minute 39)
https://stackoverflow.com/questions/20257227/how-to-find-all-shortest-paths

Input: a dictionary

    - keys label the vertices; they are assumed to be consecutive, incremental,
      nonnegative integers; the minimum vertex can be different from zero. If 
      strings are passed, they will be converted to integer before processing

    - each value is a list carrying the labels of the vertices neighboring 
      the key vertex; can be either strings or integers

'''

import numpy as np
from mpi4py import MPI
from helper_dict import changeTypeOfDictKeys
from helper_dict import renumberKeysAndValuesFrom0
from helper_json import printJson


class GeodesicEdgeBetweennessCentrality:

    def __init__(self, adjacencyList, parallel=True, sparseLabels=False):

        self.comm = MPI.COMM_WORLD
        self.me = self.comm.Get_rank()
        self.nRanks = self.comm.Get_size()
        # print("me: ", self.me, "; num procs: ", self.nRanks)

        # pre-processing of adjList 
        self.adjList = changeTypeOfDictKeys(adjacencyList,
                                            str, int)
        self.sparseLabels = sparseLabels
        self.checkForSparseLabels()

        self.nVertices = len(self.adjList)
        self.minVertex = min(list(self.adjList.keys()))
        self.gebc = np.zeros(self.nVertices, dtype=float)

        if self.me == 0:
            print(f"number of vertices: {self.nVertices}")
            print(f"min vertex: {self.minVertex}")

        if parallel:
            self.computeGebc_parallel()
        else:
            self.computeGebc_serial()


    def checkForSparseLabels(self):
        """Sometimes, a subnetwork has to be extracted from a 
        larger network, for example to obtain a network size amenable 
        to analysis. If this results in vertex ids that are no longer
        in consecutive order, some remapping becomes necessary before
        the analysis (and when results are printed)."""
        if self.sparseLabels:
            self.adjList, self.consecutiveLabelsToSparseLabels =\
                renumberKeysAndValuesFrom0(self.adjList)


    def computeGebc_serial(self):
        if self.me == 0:
            for i in range(self.minVertex, self.minVertex+self.nVertices):
                for j in range(i+1, self.minVertex+self.nVertices):
                    # this condition avoids double computations
                    self.computeGebcOfVerticesBetweenIandJ(i, j)

    def computeGebc_parallel(self):
        """use blocking communication!
        non-blocking comm leaks memory severely"""
        max = self.minVertex + self.nVertices
        if self.me == 0:
            count = 0
            n = self.nRanks - 1
            maxVertex = max-1
            for i in range(self.minVertex, max):
                if i == maxVertex:
                    # the condition i < j is no longer
                    # realized after this point
                    self.sendPoisonPill()
                    break
                destRank = count % n + 1
                self.comm.send(i, dest=destRank)
                count += 1
        else:
            i = self.comm.recv(source=0)
            while i is not None:
                for j in range(i+1, max):
                    # this condition avoids double computations
                    self.computeGebcOfVerticesBetweenIandJ(i, j)
                    # print(f"{i}, {j}", flush=True)
                i = self.comm.recv(source=0)
            print(f"rank: {self.me}, completed vertex: {i}", flush=True)
        self.aggregateGebc()

    def sendPoisonPill(self):
        for rank in range(1, self.nRanks):
            print("sending poison to: ", rank)
            self.comm.send(None, dest=rank)

    def aggregateGebc(self):
        for rank in range(1, self.nRanks):
            if self.me == 0:
                gebcFromRank = np.empty(self.nVertices, dtype=float)
                self.comm.Recv(gebcFromRank, source=rank)
                self.gebc += gebcFromRank
            else:
                self.comm.Send(self.gebc, dest=0)

    def computeGebcOfVerticesBetweenIandJ(self, i, j):
        parents, visitedParents = self.BFS(i, j)
        # print("\t---BFS done for vertices: ", i, j)
        # self.reconstructGeodesicPaths_recursive(parents, i, j, [])
        if j in parents:
            # i and j belong to
            # same connected component
            self.numPathsBetweenIandJ = 0
            self.gebcIJ = {}
            # to debug a specific pair of vertices
            # if i==5 and j==9:
            #     self.reconstructGeodesicPaths_iterative(
            #       parents, visitedParents, i, j)
            #     self.computeGebcOfVerticesInGeodesicPaths()
            self.reconstructGeodesicPaths_iterative(
                parents, visitedParents, i, j)
            self.computeGebcOfVerticesInGeodesicPaths()

    def reconstructGeodesicPaths_recursive(self, parents, startVertex,
                                           endVertex, path):
        # Parents is the list of parent vertices of each vertex,
        # as generated by BFS; proceed recursively.
        # Works fine, but highly inefficient. Not suitable for large inputs.
        # Currently broken (self.paths removed)
        path.append(endVertex)
        for parent in parents[endVertex]:
            if parent == startVertex:
                path.append(startVertex)
                self.paths.append(path.copy())
                path.pop(-1)
                break
            self.reconstructGeodesicPaths_recursive(parents, startVertex,
                                                    parent, path)
            path.pop(-1)

    def reconstructGeodesicPaths_iterative(self, parents, visited,
                                           startVertex, endVertex):
        # Parents is the list of parent vertices of each vertex
        # along a shortest route from i to j. This is equivalent
        # to a tree having root j and leaves all equal to i. Thus,
        # the problem is to traverse all branches of a tree and store
        # the corresponding path.
        # Adapted from:
        # https://www.quora.com/How-do-I-print-all-root-to-leaf-paths-in-binary-tree-iteratively
        # to debug, uncomment all print() calls

        # print("start, end: ", startVertex, endVertex)
        # print("parents: ", parents)
        # print("visited parents: ", visited)

        parents[startVertex] = []
        path = [endVertex]  # stack
        while path:
            vertex = path[-1]  # the top of the stack is the end of path
            # print("\nvertex:", vertex)
            for index, parent in enumerate(parents[vertex]):
                # print(f"parent: {parent}; index: {index}")
                if not parents[parent]:  # empty list
                    # print("Parent is a leaf (startVertex)."
                    #       + "\n\t************ Print path:", path)
                    # print("\tPop ", vertex)
                    self.updateGebcOfVerticesInGeodesicPath(path, endVertex)
                    visited[vertex][index] = True
                    path.pop()
                    break
                if not visited[vertex][index]:
                    path.append(parent)
                    visited[vertex][index] = True
                    # print(f"{parent} pushed on the path")
                    break
                if parent == parents[vertex][-1]:
                    path.pop()
                    visited[vertex] = [
                        False for _ in range(len(visited[vertex]))
                    ]
                    # print(f"All parents of vertex visited."
                    #       + f"\n\tPop {vertex} from path.")
                    break
            # print("visited:", visited)

    def updateGebcOfVerticesInGeodesicPath(self, path, endVertex):
        self.numPathsBetweenIandJ += 1
        for vertex in path:
            if vertex not in self.gebcIJ:
                self.gebcIJ[vertex] = 1
            else:
                self.gebcIJ[vertex] += 1
        # Start and end vertices trivially belong to any
        # geodesic path connecting them, so do not count
        # these contributions. StartVertex was never added
        # to path, while endVertex was (always exactly once)
        self.gebcIJ[endVertex] -= 1

    def computeGebcOfVerticesInGeodesicPaths(self):
        # add the contribution of i and j to gebc of vertex
        for vertex in self.gebcIJ.keys():
            self.gebc[vertex] += self.gebcIJ[vertex]/self.numPathsBetweenIandJ

    # breadth first search
    def BFS(self, startVertex, endVertex):

        # Find all shortest paths from startVertex to endVertex.
        # At this stage, encode this information in a list containing
        # the parent vertices of each vertex on each path.

        parents = {}
        visitedParents = {}  # for geodesic path reconstruction later
        distance = {}  # from startVertex
        queue = []
        visited = np.zeros(self.nVertices, dtype=int)

        queue.append(startVertex)
        visited[startVertex] = 1
        distance[startVertex] = 0
        parents[startVertex] = [startVertex]
        visitedParents[startVertex] = [False]

        while queue:
            currentVertex = queue.pop(0)
            if currentVertex == endVertex:
                break
            for adjVertex in self.adjList[currentVertex]:
                if not visited[adjVertex]:
                    queue.append(adjVertex)
                    distance[adjVertex] = distance[currentVertex] + 1
                    parents[adjVertex] = [currentVertex]
                    visitedParents[adjVertex] = [False]
                    visited[adjVertex] = 1
                # If the neighbor adjVertex has been visited already,
                # then it has a parent. If this parent has the same distance
                # from startVertex as currentVertex, then currentVertex is also
                # a parent of adjVertex
                elif distance[parents[adjVertex][0]] == distance[currentVertex]:
                    parents[adjVertex].append(currentVertex)
                    visitedParents[adjVertex].append(False)
        return parents, visitedParents

    def printJson(self, fileName, normalize=True):
        """
        If sparse vertex ids were remapped to consecutive 
        vertex ids, apply the back-conversion map to the 
        indices of self.gebc, then create a dict and print 
        it in json format.
        """
        if self.me == 0:
            if normalize:
                self.gebc /= max(self.gebc)
            gebc = {}
            if self.consecutiveLabelsToSparseLabels:
                for index, value in enumerate(self.gebc):
                    gebc[self.consecutiveLabelsToSparseLabels[index]] = value
            else:
                for index, value in enumerate(self.gebc):
                    gebc[index] = value
            printJson(gebc, fileName)


    def closeMPI(self):
        MPI.Finalize()
