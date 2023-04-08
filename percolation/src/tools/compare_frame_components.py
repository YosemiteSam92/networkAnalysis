#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))
from graph.graph_structs import MetaGraph, EdgeTag
from mols.trajectorySnapshot import TrajectorySnapshot
import numpy as np
import argparse
import resource


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
        "-c1",
        "--connectivity1",
        default="data/first_conn.dat",
        help="The path to the analysis connectivity file.",
    )
    parser.add_argument(
        "-t1",
        "--trajectory1",
        default="data/first_traj.dat",
        help="The path to the analysis trajectory file.",
    )

    parser.add_argument(
        "-c2",
        "--connectivity2",
        default="data/second_conn.dat",
        help="The path to the analysis connectivity file.",
    )
    parser.add_argument(
        "-t2",
        "--trajectory2",
        default="data/second_conn.traj",
        help="The path to the analysis trajectory file.",
    )

    args = parser.parse_args()

    # tena2 paths
    dumpConnectivity1 = args.connectivity1
    dumpTrajectory1 = args.trajectory1
    dumpConnectivity2 = args.connectivity2
    dumpTrajectory2 = args.trajectory2

    # tena2 only
    resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))
    sys.setrecursionlimit(10 ** 6)

    data = None

    connectivity_reader = SingleRecordReader(data)
    trajectory_reader = SingleRecordReader(data)
    aggregate_reader = AggregateRecordReader(
        data, [connectivity_reader, trajectory_reader]
    )

    graph1 = None
    graph2 = None

    components1 = None
    components2 = None

    with open(dumpConnectivity1, "r") as connectivity_file:
        with open(dumpTrajectory1, "r") as trajectory_file:
            aggr_record = aggregate_reader.readRecord(
                [connectivity_file, trajectory_file]
            )

            graph1 = build_graph(aggr_record)

            num_comps, comp_data = graph1.get_components()

            components1 = []

            for comp in comp_data:
                comp_nodes = []
                for a, d in comp:
                    comp_nodes.append(a)
                components1.append(set(comp_nodes))

    with open(dumpConnectivity2, "r") as connectivity_file:
        with open(dumpTrajectory2, "r") as trajectory_file:
            aggr_record = aggregate_reader.readRecord(
                [connectivity_file, trajectory_file]
            )

            graph2 = build_graph(aggr_record)

            num_comps, comp_data = graph2.get_components()

            components2 = []

            for comp in comp_data:
                comp_nodes = []
                for a, d in comp:
                    comp_nodes.append(a)
                components2.append(set(comp_nodes))

    for seta in components1:
        for setb in components2:
            intersection = seta.intersection(setb)

            if intersection and not seta.issubset(setb):
                print("Elements of the following subset are no longer connected:")
                print(seta)
                print("The following elements have lost connection to the rest")
                print(seta - setb)
