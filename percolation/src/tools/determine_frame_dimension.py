#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))
from graph.graph_structs import MetaGraph, EdgeTag
from mols.trajectorySnapshot import TrajectorySnapshot
import numpy as np


def splitBondEntry(l):
    line = l.split(" ")
    bondtype, partner_1, partner_2 = int(line[1]), int(line[2]), int(line[3])
    return bondtype, partner_1, partner_2


def build_graph(aggr_record):
    bondRecord, trajRecord = aggr_record.record_data
    print("timestep:", trajRecord.timeStep)
    traj = TrajectorySnapshot(trajRecord)
    
    graph = MetaGraph()
    
    xpos = traj.x
    
    for bond_entry in bondRecord.entries:
        bondtype, a1, a2 = splitBondEntry(bond_entry)
        
        dist = xpos[a2]-xpos[a1]
        
        raw_tag = []
        
        for i in range(len(dist)):
            if dist[i] < -0.5: 
                raw_tag.append(-1)
            elif dist[i] > 0.5:
                raw_tag.append(1)
            else:
                raw_tag.append(0)
        
        tag = EdgeTag(*raw_tag)
        
        graph.add_edge(a1,a2, tag)    

    return graph


import argparse
import resource

if __name__ == "__main__":
    from streaming.LammpsData import LammpsData
    from streaming.hashForTimeOrdering import HashForTimeOrdering
    from mols.molecularGraph import MolecularGraph
    from mols.trajectorySnapshot import TrajectorySnapshot
    from streaming.dump_reader import (
        SingleRecordReader,
        AggregateRecordReader,
        MultipleRecordReader,
    )
    from workers.recordProcessor import MultipleRecordProcessor

    parser = argparse.ArgumentParser(
        description="Script to analyze dimension of one frame of a trajectory and connectivity set."
    )
    parser.add_argument(
        "-c",
        "--connectivity",
        default="/home/mattia/Documents/emmy_backup/large_503K_3/dump/dump.1_btype_batom1_batom2_nptCuring_gelation",
        help="The path to the analysis connectivity file.",
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        default="/home/mattia/Documents/emmy_backup/large_503K_3/dump/dump.1_atomId_type_scaledCoords_nptCuring_gelation",
        help="The path to the analysis trajectory file.",
    )

    args = parser.parse_args()

    # tena2 paths
    dumpConnectivity = args.connectivity
    dumpTrajectory = args.trajectory
    
    # tena2 only
    resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))
    sys.setrecursionlimit(10 ** 6)
    
    data = None

    connectivity_reader = SingleRecordReader(data)
    trajectory_reader = SingleRecordReader(data)
    aggregate_reader = AggregateRecordReader(
        data, [connectivity_reader, trajectory_reader]
    )

    with open(dumpConnectivity, "r") as connectivity_file:
        with open(dumpTrajectory, "r") as trajectory_file:
            aggr_record = aggregate_reader.readRecord([connectivity_file, trajectory_file])
            
            graph = build_graph(aggr_record)
            
            num_comps, dims = graph.find_stable_loops()
            
            max_dim = np.max(np.array(dims))

            print("Frame has maximum grid dimension {0}".format(max_dim))
