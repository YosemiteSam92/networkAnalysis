import argparse
import resource
import sys

if __name__ == "__main__":
    import sys
    from os.path import dirname, abspath

    sys.path.append(abspath(dirname(__file__)))

    from streaming.LammpsData import LammpsData
    from mols.molecularGraph import MolecularGraph
    from mols.largestMolecule import LargestMolecule
    from mols.trajectorySnapshot import TrajectorySnapshot
    from streaming.dump_reader import (
        SingleRecordReader,
        AggregateRecordReader,
        MultipleRecordReader,
    )
    from workers.recordProcessor import MultipleRecordProcessor

    parser = argparse.ArgumentParser(
        description="Script to print the atomic coordinates of the largest"
                    + "molecule and of all remaining molecules in two separate xyz files, "
                    + "for each time frame."
    )
    parser.add_argument(
        "-d",
        "--datafile",
        default="/home/mattia/Documents/emmy_backup/large_503K_1/data/data.Min",
        help="path to lammps data file",
    )
    parser.add_argument(
        "-c",
        "--connectivity",
        default="/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump.1_btype_batom1_batom2_nptCuring_gelation",
        help="path to connectivity dump file",
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        default="/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump.1_atomId_type_scaledCoords_nptCuring_gelation",
        help="path to trajectory dump file",
    )
    parser.add_argument(
        "-oL",
        "--outputLargest",
        default="/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large1/curing/largestMolecules",
        help="directory that will host a single xyz file for the largest molecule, for each frame",
    )
    parser.add_argument(
        "-oO",
        "--outputOthers",
        default="/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large1/curing/allOtherMolecules",
        help="directory that will host a single xyz file for all but the largest molecules, for each frame",
    )

    args = parser.parse_args()


    def printLargestMoleculeAndEverythingElseToXyz(aggr_record, varargs):
        bondRecord, trajRecord = aggr_record
        print("timestep:", trajRecord.timeStep)
        traj = TrajectorySnapshot(trajRecord)
        molGraph = MolecularGraph(bondRecord, traj)
        molGraph.computeMolecularMasses_numOfIntramolCrosslinks_largestMolecule()
        largestMolecule = LargestMolecule(molGraph)
        outputFiles = varargs[1]
        largestMolecule.printToXyz(outputFiles[0], outputFiles[1])


    dataFile = args.datafile
    dumpConnectivity = args.connectivity
    dumpTrajectory = args.trajectory
    outputFiles = (args.outputLargest, args.outputOthers)

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

            # This is the intended way to pass them, I just screwed it up
            multi_processor.processEntriesInParallel(
                printLargestMoleculeAndEverythingElseToXyz,
                outputFiles 
                # put here any further argument to the record-
                # processing function 
            )
