#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))
from mols.molecularGraph import MolecularGraph
from mols.trajectorySnapshot import TrajectorySnapshot


import argparse
import resource


def splitBondEntry(l):
    line = l.split(" ")
    bondtype, partner_1, partner_2 = int(line[1]), int(line[2]), int(line[3])
    return bondtype, partner_1, partner_2


if __name__ == "__main__":
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
        "-c1",
        "--connectivity1",
        help="The path to the first analysis connectivity file.",
    )
    parser.add_argument(
        "-c2",
        "--connectivity2",
        help="The path to the second analysis connectivity file.",
    )

    args = parser.parse_args()

    # tena2 paths
    dumpConnectivity1 = args.connectivity1
    dumpConnectivity2 = args.connectivity2

    # tena2 only
    resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))
    sys.setrecursionlimit(10 ** 6)

    connectivity_reader = SingleRecordReader(None)

    first_bonds = []
    first_pairs = []
    second_bonds = []
    second_pairs = []

    with open(dumpConnectivity1, "r") as connectivity_file:
        record = connectivity_reader.readRecord(connectivity_file)

        for entry in record.entries:
            bondtype, partner_1, partner_2 = splitBondEntry(entry)
            first_bonds.append((bondtype, partner_1, partner_2))
            first_bonds.append((bondtype, partner_2, partner_1))
            first_pairs.append((partner_1, partner_2))
            first_pairs.append((partner_2, partner_1))

    with open(dumpConnectivity2, "r") as connectivity_file:
        record = connectivity_reader.readRecord(connectivity_file)

        for entry in record.entries:
            bondtype, partner_1, partner_2 = splitBondEntry(entry)
            second_bonds.append((bondtype, partner_1, partner_2))
            second_bonds.append((bondtype, partner_2, partner_1))
            second_pairs.append((partner_1, partner_2))
            second_pairs.append((partner_2, partner_1))

    fb_set = set(first_bonds)
    fp_set = set(first_pairs)

    sb_set = set(second_bonds)
    sp_set = set(second_pairs)

    first_not_second_bonds = fb_set.difference(sb_set)
    first_not_second_pair = fp_set.difference(sp_set)

    second_not_first_bonds = sb_set.difference(fb_set)
    second_not_first_pair = sp_set.difference(fp_set)

    print(len(first_not_second_bonds), len(second_not_first_bonds))
    print("First, not second:")
    for t, a, b in first_not_second_bonds:
        print("{0} --> {1} [{2}]".format(a, b, t))

    print("")
    print("Second, not first:")
    for t, a, b in second_not_first_bonds:
        print("{0} --> {1} [{2}]".format(a, b, t))

    print(len(first_not_second_pair), len(second_not_first_pair))
    print("First, not second:")
    for a, b in first_not_second_pair:
        print("{0} --> {1}".format(a, b))

    print("")
    print("Second, not first:")
    for a, b in second_not_first_pair:
        print("{0} --> {1} ".format(a, b))
