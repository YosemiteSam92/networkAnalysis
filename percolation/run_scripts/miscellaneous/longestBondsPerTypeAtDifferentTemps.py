#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "../"))
from streaming.singleRecord_indexBtypeBlength import SingleRecord_indexBtypeBlength

raise Exception("Well, shit, this is still broken :)\n\tMy bad ^^")

from anlammps.multipleDumpsWithSingleRecord import MultipleDumpsWithSingleRecord

def printProcessedRecords(multipleDumpsWithSingleRecord, outputFiles):
    i = 0
    for processedSingleRecord in multipleDumpsWithSingleRecord.processedSingleRecords:
        outputFile = outputFiles[i]
        processedSingleRecord.printBondLengthsOfDifferentBondTypesInDIfferentColumns(
            outputFile
        )
        i += 1


dumpFiles = [
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/dump.300C_btype_dist",
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/dump.350C_btype_dist",
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/dump.400C_btype_dist",
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/dump.450C_btype_dist",
]

outputFiles = [
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/300C_blengthsInColumnsByBtype",
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/350C_blengthsInColumnsByBtype",
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/400C_blengthsInColumnsByBtype",
    "/Users/mattia/Documents/epoxy/MDsystems/large2/dumpBtypeDist/450C_blengthsInColumnsByBtype",
]

multipleDumpsWithSingleRecord = MultipleDumpsWithSingleRecord(dumpFiles, None)
multipleDumpsWithSingleRecord.processSingleRecords(SingleRecord_indexBtypeBlength)
printProcessedRecords(multipleDumpsWithSingleRecord, outputFiles)
