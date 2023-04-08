#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))
from graph.graph_structs import MetaGraph, EdgeTag
from mols.trajectorySnapshot import TrajectorySnapshot
import numpy as np
import argparse

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

        dist = xpos[a2] - xpos[a1]

        raw_tag = []

        for i in range(len(dist)):
            if dist[i] < -0.5:
                raw_tag.append(-1)
            elif dist[i] > 0.5:
                raw_tag.append(1)
            else:
                raw_tag.append(0)

        tag = EdgeTag(*raw_tag)

        graph.add_edge(a1, a2, tag)

    return graph

def get_graph_dim(graph):
    num_comps, dims = graph.find_stable_loops()
    
    max_dim = np.max(np.array(dims))
    return max_dim

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
        description="Script to dump the meta graph representation of frame data."
    )
    parser.add_argument(
        "-c",
        "--connectivity",
        default="data/first_conn.dat",
        help="The path to the analysis connectivity file.",
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        default="data/first_traj.dat",
        help="The path to the analysis trajectory file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="data/frame_",
        help="The path prefix to the output dump files.",
    )

    args = parser.parse_args()

    dumpConnectivity1 = args.connectivity
    dumpTrajectory1 = args.trajectory
    outputPath = args.output

    data = None

    connectivity_reader = SingleRecordReader(data)
    trajectory_reader = SingleRecordReader(data)
    aggregate_reader = AggregateRecordReader(
        data, [connectivity_reader, trajectory_reader]
    )

    graph1 = None

    with open(dumpConnectivity1, "r") as connectivity_file:
        with open(dumpTrajectory1, "r") as trajectory_file:
            aggr_record = aggregate_reader.readRecord(
                [connectivity_file, trajectory_file]
            )

            graph1 = build_graph(aggr_record)

    graph1.mark_states()
    graph1.dump(outputPath + "_full.dot", True)
    print("Raw:",get_graph_dim(graph1))

    """clone1 = graph1.copy()
    clone1 = clone1.get_simplified_graph()
    clone1 = clone1.get_purged_graph()
    
    print("Simplified purged:",get_graph_dim(clone1))
    

    clone1 = graph1.copy()
    clone1.reduce()
    clone1 = clone1.get_purged_graph()
    
    print("Reduced:",get_graph_dim(clone1))"""
    
    graph1.reduce()
    simple1 = graph1.get_simplified_graph()
    simple1.dump(outputPath + "_simple.dot", True)
    print("Simplified:",get_graph_dim(simple1))
    
    purged1 = simple1.get_purged_graph()
    purged1.dump(outputPath + "_simple_purge.dot", True)    
    print("Simplified purged:",get_graph_dim(purged1))
    
    comp_graph = purged1.get_component_graph()
    print("Component:",get_graph_dim(comp_graph))
    
    comp_graph.reduce()
    comp_graph.dump(outputPath + "_components.dot", True)
    print("Component reduced:",get_graph_dim(comp_graph))

    purged = comp_graph.get_purged_graph()
    purged.mark_states()
    purged.dump(outputPath + "_components_purged.dot", True)
    print("Component reduced purged:",get_graph_dim(purged))

    purged_simplified = purged.get_simplified_graph()
    purged_simplified.mark_states()
    purged_simplified.dump(outputPath + "_components_purged_simplified.dot", True)
    print("Component reduced purged simplified:",get_graph_dim(purged_simplified))
    
    purged_simplified2 = purged_simplified.get_simplified_graph()
    purged_simplified2.mark_states()
    purged_simplified2.dump(outputPath + "_components_purged_simplified_2.dot", True)
    print("Component reduced purged simplified x 2:",get_graph_dim(purged_simplified2))