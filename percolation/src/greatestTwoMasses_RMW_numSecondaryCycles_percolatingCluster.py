#!/usr/bin/env python3

import argparse
import resource

if __name__ == "__main__":
    import sys
    from os.path import dirname, abspath

    sys.path.append(abspath(dirname(__file__)))
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
        default="/home/mattia/Documents/emmy_backup/large_503K_1/data/data.Min",
        help="The path to the analysis data file.",
    )
    parser.add_argument(
        "-c",
        "--connectivity",
        default="/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump.1_btype_batom1_batom2_nptCuring_gelation",
        help="The path to the analysis connectivity file.",
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        default="/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump.1_atomId_type_scaledCoords_nptCuring_gelation",
        help="The path to the analysis trajectory file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large1/gelPoint/gelPointMeasurements1_2nsNptCuring",
        help="Path to output file.",
    )

    args = parser.parse_args()

    def computeGreatestMolecularMasses_RMW_numOfIntramolCrosslinks_existenceOfPercolatingCluster(
        aggr_record, varargs
    ):
        bondRecord, trajRecord = aggr_record
        print("timestep:", trajRecord.timeStep)
        traj = TrajectorySnapshot(trajRecord)
        molGraph = MolecularGraph(bondRecord, traj)
        molGraph.computeMolecularMasses_numOfIntramolCrosslinks_largestMolecule()
        molGraph.computeRMW()
        molGraph.getTwoGreatestMolecularMasses()
        molGraph.isAnyMoleculePercolating()
        outputQueueEntry = (
            molGraph.greatestTwoMasses[0],
            molGraph.greatestTwoMasses[1],
            molGraph.numIntramolCrosslinks,
            molGraph.RMW,
            molGraph.numChains,
            molGraph.numSheets,
            molGraph.numPercolatingMolecules,
            molGraph.largestMoleculePercolating,
        )
        return outputQueueEntry

    def printOutputQueue(outputQueue, outputFile, outputHeader):
        print("Printing output queue")
        hash = HashForTimeOrdering(outputQueue)
        hash.printToFile(outputFile, outputHeader)

    # dataFile = "/Users/mattia/Documents/anlammps/lammps/lammpsData/data.Min"
    # dumpConnectivity = "/Users/mattia/Documents/anlammps/lammps/lammpsDump/dump.1_btype_batom1_batom2_nptCuring_gelation_short"
    # dumpTrajectory = "/Users/mattia/Documents/anlammps/lammps/lammpsDump/dump.1_atomId_type_scaledCoords_nptCuring_gelation_short"
    # outputFile = "/Users/mattia/Documents/anlammps/scripts/testOutput/RMW_secondaryCycles_percolatingCluster.txt"

    # tena2 paths
    dataFile = args.datafile
    dumpConnectivity = args.connectivity
    dumpTrajectory = args.trajectory
    outputFile = args.output

    # tena2 only
    resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))
    sys.setrecursionlimit(10 ** 6)

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

            multi_processor = MultipleRecordProcessor(sequential_reader, data)

            multi_processor.processEntriesInParallel(
                computeGreatestMolecularMasses_RMW_numOfIntramolCrosslinks_existenceOfPercolatingCluster,
                [],
            )

            outputHeader = (
                "timestep greatestMass secondGreatestMass numOfIntramolCrosslinks"
                + " RMW numChains numSheets numOfPercolatingMolecules isLargestMoleculePercolating"
            )
            printOutputQueue(multi_processor.outputQueue, outputFile, outputHeader)
