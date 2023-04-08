import argparse
import resource
import sys

if __name__ == "__main__":
    import sys
    from os.path import dirname, abspath

    sys.path.append(abspath(dirname(__file__)))

    from streaming.LammpsData import LammpsData
    from streaming.hashForTimeOrdering import HashForTimeOrdering
    from mols.molecularGraph import MolecularGraph
    from streaming.dump_reader import (
        SingleRecordReader,
        MultipleRecordReader,
    )
    from workers.recordProcessor import MultipleRecordProcessor

    parser = argparse.ArgumentParser(
        description="Script to analyze the greatest two masses"
        + "the reduced molecular weight and the number"
        + "of secondary cycles."
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
        "-o",
        "--output",
        default="/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large1/gelPoint/gelPointMeasurements1_2nsNptCuring",
        help="Path to output directory.",
    )

    args = parser.parse_args()

    def computeGreatestMolecularMasses_RMW_numOfIntramolCrosslinks(bondRecord, varargs):
        print("timestep:", bondRecord.get_timestep())
        molGraph = MolecularGraph(bondRecord, [])
        molGraph.computeMolecularMasses_numOfIntramolCrosslinks_largestMolecule()
        molGraph.computeRMW()
        molGraph.getTwoGreatestMolecularMasses()
        outputQueueEntry = (
            molGraph.greatestTwoMasses[0],
            molGraph.greatestTwoMasses[1],
            molGraph.numIntramolCrosslinks,
            molGraph.RMW,
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
    outputFile = args.output

    # tena2 only
    resource.setrlimit(resource.RLIMIT_STACK, (2 ** 29, -1))
    sys.setrecursionlimit(10 ** 6)

    data = LammpsData(dataFile, ["atoms"])

    connectivity_reader = SingleRecordReader(data)
    sequential_reader = MultipleRecordReader(data, connectivity_reader)

    with open(dumpConnectivity, "r") as connectivity_file:
        sequential_reader.setDumpfile(connectivity_file)

        multi_processor = MultipleRecordProcessor(sequential_reader, data)

        multi_processor.processEntriesInParallel(
            computeGreatestMolecularMasses_RMW_numOfIntramolCrosslinks,
            [],
        )

        outputHeader = (
            "timestep greatestMass secondGreatestMass numOfIntramolCrosslinks" + " RMW"
        )
        printOutputQueue(multi_processor.outputQueue, outputFile, outputHeader)
