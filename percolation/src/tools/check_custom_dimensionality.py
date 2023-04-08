#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))
import argparse
from graph.graph_structs import MetaGraph, EdgeTag
from mols.trajectorySnapshot import TrajectorySnapshot
import numpy as np


def splitBondEntry(l):
    line = l.split(" ")
    bondtype, partner_1, partner_2 = int(line[1]), int(line[2]), int(line[3])
    return bondtype, partner_1, partner_2


def get_graph_dim(aggr_record):
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
        description="Script to analyze the msd and diffusion coefficients of buckets."
    )
    parser.add_argument(
        "-d",
        "--datafile",
        default="/home/mattia/Documents/scriptDeveloper/anlammps/lammps/lammpsData/data.Min_large3",
        help="The path to the analysis data file.",
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
    parser.add_argument(
        "-o",
        "--output",
        default="./data/check_",
        help="Output path prefix.",
    )

    args = parser.parse_args()

    # tena2 paths
    dataFile = args.datafile
    dumpConnectivity = args.connectivity
    dumpTrajectory = args.trajectory
    output_prefix = args.output

    data = LammpsData(dataFile, ["atoms"])

    connectivity_reader = SingleRecordReader(data)
    trajectory_reader = SingleRecordReader(data)
    aggregate_reader = AggregateRecordReader(
        data, [connectivity_reader, trajectory_reader]
    )
    sequential_reader = MultipleRecordReader(data, aggregate_reader)

    with open(dumpConnectivity, "r") as connectivity_file:
        with open(dumpTrajectory, "r") as trajectory_file:
            sequential_reader.setDumpfile([connectivity_file, trajectory_file])

            curr_dim = -1
            last_record = None

            while True:
                next_record = sequential_reader.readRecord()
                if next_record is None:
                    break
                new_dim = get_graph_dim(next_record)

                if new_dim < curr_dim:
                    last_record.dump_to(
                        [
                            output_prefix + "first_conn.dat",
                            output_prefix + "first_traj.dat",
                        ]
                    )
                    next_record.dump_to(
                        [
                            output_prefix + "second_conn.dat",
                            output_prefix + "second_traj.dat",
                        ]
                    )
                    print(
                        "Dimension {0} exceeds succeeding dimension {1}".format(
                            curr_dim, new_dim
                        )
                    )
                    break

                if new_dim != curr_dim:
                    print("Now at dimension {0}".format(new_dim))
                    curr_dim = new_dim
                last_record = next_record
