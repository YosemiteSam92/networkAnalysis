from multiprocessing import Queue


def sortKeysInAscendingOrder(dict):
    keys = list(dict.keys())
    keys.sort()
    return keys


def printHeader(output, header):
    print(header, file=output)


def joinTupleOfFloatsWithBlanks(tuple_data):
    result_string = " ".join([str(x) for x in tuple_data])
    return result_string


# Class to print data from a queue ordered by time
class HashForTimeOrdering:
    def __init__(self, outputQueue: Queue):
        # outputQueue is a multiprocessing.Queue object
        # its first entry is always record.timeStep
        # its second entry is a tuple, containing various output quantities
        self.initAttributes(outputQueue)
        self.populateHash()
        self.sortedTimeKeys = sortKeysInAscendingOrder(self.hash)

    def initAttributes(self, outputQueue):
        self.hash = {}
        self.outputQueue = outputQueue

    def populateHash(self):
        while 1:
            try:
                self.createHashEntry()
            except:
                break

    def createHashEntry(self):
        key, outputTuple = self.outputQueue.get(False)
        # print("from queue:", key)
        self.hash[key] = outputTuple

    def printToFile(self, outputFile, header):
        with open(outputFile, "w") as output:
            printHeader(output, header)
            self.printAllSortedKeyTuplePairs(output)

    def printAllSortedKeyTuplePairs(self, output):
        for timeKey in self.sortedTimeKeys:
            self.printThisKeyTuplePair(timeKey, output)

    def printThisKeyTuplePair(self, timeKey, output):
        line = self.createOutputLine(timeKey)
        print(line, file=output)

    def createOutputLine(self, timeKey):
        line = format(timeKey) + " "  # convert to string
        line += joinTupleOfFloatsWithBlanks(self.hash[timeKey])
        return line
