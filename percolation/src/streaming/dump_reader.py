#!/usr/bin/env python3


# Class holding the data of a single record
class SingleRecord:
    def __init__(self, data):
        # Additional data for the entry (opaque)
        self.data = data
        # For dumping: Line containing item format information
        self.item_format_line = None
        # Time step identifier
        self.timeStep = None
        # Number of entries
        self.numEntries = None
        # Bounding box information
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.zhi = None

        # actual entries in record
        self.entries = []

    def get_timestep(self):
        return self.timeStep

    def print_info(self):
        print("Info of Single Record:")
        print("Time Step:", self.timeStep)
        print("Num entries:", self.numEntries)
        print("Bounds:")
        print(self.xlo, "\t", self.xhi)
        print(self.ylo, "\t", self.yhi)
        print(self.zlo, "\t", self.zhi)
        print("Format:")
        print(self.item_format_line)
        print("Actual number of entries:", len(self.entries))

    def dump_to(self, filename):
        with open(filename, "w") as out:
            out.write("ITEM: TIMESTEP\n")
            out.write("{0}\n".format(self.timeStep))
            out.write("ITEM: NUMBER OF ATOMS\n")
            out.write("{0}\n".format(self.numEntries))
            out.write("ITEM: BOX BOUNDS pp pp pp\n")
            out.write("{0} {1}\n".format(self.xlo, self.xhi))
            out.write("{0} {1}\n".format(self.ylo, self.yhi))
            out.write("{0} {1}\n".format(self.zlo, self.zhi))
            out.write("{0}\n".format(self.item_format_line))
            for rec in self.entries:
                out.write("{0}\n".format(rec))


class AggregateRecord:
    def __init__(self, record_data, data):
        self.data = data
        self.record_data = record_data

    def __getitem__(self, key):
        if key < 0 or key > len(self.record_data):
            return None

        return self.record_data[key]

    def get_timestep(self):
        if self.record_data:
            return self.record_data[0].get_timestep()
        else:
            return -1

    def print_info(self):
        print("Info of Aggregate Record:")
        print("[")
        for e in self.record_data:
            e.print_info()
            print(",")
        print("]")

    def dump_to(self, filenames):
        if len(filenames) != len(self.record_data):
            raise Exception("Wrong number of file paths to dump to")

        for i in range(len(filenames)):
            self.record_data[i].dump_to(filenames[i])


def splitLine(line):
    return line.split(" ")


# warning: round does not always round the last digit correctly
def readLo(bounds):
    return round(float(bounds[0]), 3)


def readHi(bounds):
    return round(float(bounds[1]), 3)


def readNonEmptyLine(file):
    while 1:
        line = file.readline()
        # deal with end of file
        if line == "":
            return None

        if line.strip() != "":
            return line.rstrip(" \r\n")


# input: file handle pointing to the first line of a dump file, with a single
# or multiple records
# output: self.record, a list containing all of the record entries included
# between "ITEM: ENTRIES ..." (excluded) and the end of the record; each line is
# a separate list entry, stored as a string deprived of any terminating blank
# and/or EOL characters.
#
# Two possible formats of the entries:
#
# entry format: index bondType bondAtom1 bondAtom2
#
# entry format: atomid atomType x y z
# coordinates are scaled
class SingleRecordReader:
    def __init__(self, data):
        self.data = data
        self.error = False

    def readRecord(self, dump):
        result = SingleRecord(self.data)
        # Read Preambol

        # Skip timestep ITEM line
        dump.readline()

        # Read time step
        timeStepLine = readNonEmptyLine(dump)
        if timeStepLine is None:
            self.error = True
            return None
        result.timeStep = int(timeStepLine)

        # Skip NumEntries ITEM line
        dump.readline()

        # Read number of entries
        numEntriesLine = readNonEmptyLine(dump)
        if numEntriesLine is None:
            self.error = True
            return None
        result.numEntries = int(numEntriesLine)

        # Skip Bounding box ITEM line
        dump.readline()

        # Read bounds of box
        boundsLine = readNonEmptyLine(dump)
        if boundsLine is None:
            self.error = True
            return None

        bounds = splitLine(boundsLine)
        result.xlo = readLo(bounds)
        result.xhi = readHi(bounds)

        boundsLine = readNonEmptyLine(dump)
        if boundsLine is None:
            result.error = True
            return None

        bounds = splitLine(boundsLine)
        result.ylo = readLo(bounds)
        result.yhi = readHi(bounds)

        boundsLine = readNonEmptyLine(dump)
        if boundsLine is None:
            self.error = True
            return None

        bounds = splitLine(boundsLine)
        result.zlo = readLo(bounds)
        result.zhi = readHi(bounds)

        # Read line containing entry format description
        result.item_format_line = readNonEmptyLine(dump)

        # Read entries as record
        for i in range(result.numEntries):
            line = readNonEmptyLine(dump)

            if line is None:
                self.error = True
                return None

            result.entries.append(line)

        return result


# Class to aggregate records from multiple input streams into one
class AggregateRecordReader:
    def __init__(self, data, entryReaders=SingleRecordReader):
        self.data = data
        self.error = False

        if type(entryReaders) is list:
            self.entryReaders = entryReaders
        else:
            self.entryReaders = [entryReaders]

    def readRecord(self, dumpFiles):
        if type(dumpFiles) is list:
            dumpFiles = dumpFiles
        else:
            dumpFiles = [dumpFiles]

        if len(self.entryReaders) != len(dumpFiles):
            self.error = True
            raise Exception("incompatible number of files and readers")

        total_record = []
        for i in range(len(self.entryReaders)):
            sub_entry = self.entryReaders[i].readRecord(dumpFiles[i])
            if sub_entry is None or self.entryReaders[i].error:
                self.error = True
                return None
            # sub_entry.print_info()
            total_record.append(sub_entry)

        # Return the entry
        return AggregateRecord(total_record, self.data)


# Class to read multiple records from one file by emplyoing a reader for individual records
class MultipleRecordReader:
    # a record contains all of the info
    # relative to a single timestep; it is made of multiple entries

    def __init__(self, data, entryReader=SingleRecordReader, num_entries=0):
        self.total_entries = num_entries
        self._entries_read = 0
        self._entryReader = entryReader
        self.error = False
        self.dumpFile = None

    def setDumpfile(self, dumpFile):
        self.dumpFile = dumpFile

    def readRecord(self, dumpFile=None):
        if self.total_entries > 0 and self._entries_read >= self.total_entries:
            # Do not read past the last entry
            return None
        else:
            if dumpFile is None:
                entry = self._entryReader.readRecord(self.dumpFile)
            else:
                entry = self._entryReader.readRecord(dumpFile)

            if entry is None or self._entryReader.error:
                self.error = True
                return None

            self._entries_read += 1

            # Return the entry
            return entry
