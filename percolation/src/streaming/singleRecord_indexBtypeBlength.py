import sys


class SingleRecord_indexBtypeBlength:
    # expect record of format: index Btype Blength
    def __init__(self, record):
        self.record = record
        self.bondLenghtsByBondType = {}
        # key: bond type
        # value: a list containing all lengths of bonds with type "bond type"
        self.sortBondLengthsByBondType()

    def sortBondLengthsByBondType(self):
        for entry in self.record:
            btype, blength = splitEntry(entry)
            self.populateBondLenghtsByBondType(btype, blength)

    def populateBondLenghtsByBondType(self, btype, blength):
        btypeInt = int(btype)
        if btypeInt in self.bondLenghtsByBondType:
            self.bondLenghtsByBondType[btypeInt].append(float(blength))
        else:
            self.bondLenghtsByBondType[btypeInt] = [float(blength)]

    ######################      print methods    ######################

    def printBondLengthsOfDifferentBondTypesInDIfferentColumns(self, outputFile):
        with open(outputFile, "w") as output:
            self.printHeader(output)
            self.printLinesWithBondLengthsOfDifferentTypesSeparatedByBlanks(output)

    def printHeader(self, output):
        print(self.createHeader(), file=output)

    def createHeader(self):
        self.sortBondTypesAsInt()
        return "\t".join(map(str, self.sortedBondTypesAsInt))

    def sortBondTypesAsInt(self):
        self.sortedBondTypesAsInt = list(self.bondLenghtsByBondType.keys())
        self.sortedBondTypesAsInt.sort()

    def printLinesWithBondLengthsOfDifferentTypesSeparatedByBlanks(self, output):
        self.findBondTypeWithGreatestNumberOfBonds()
        for rowNum in range(self.greatestNumOfBondsOfACertainType):
            self.printLineWithBondLengthsOfDifferentTypesSeparatedByBlanks(rowNum, output)

    def findBondTypeWithGreatestNumberOfBonds(self):
        numBondsOfEachTypeSortedByIncreasingBondType = self.storeNumOfBondsOfEachType()
        self.greatestNumOfBondsOfACertainType = max(numBondsOfEachTypeSortedByIncreasingBondType)

    def storeNumOfBondsOfEachType(self):
        numOfBondsOfEachTypeSortedByIncreasingBondType = []
        for btypeInt in self.sortedBondTypesAsInt:
            numOfBondsOfEachTypeSortedByIncreasingBondType.append(len(self.bondLenghtsByBondType[btypeInt]))
        return numOfBondsOfEachTypeSortedByIncreasingBondType

    def printLineWithBondLengthsOfDifferentTypesSeparatedByBlanks(self, rowNum, output):
        print(self.lineWithBondLengthsOfDifferentTypesSeparatedByBlanks(rowNum), file=output)

    def lineWithBondLengthsOfDifferentTypesSeparatedByBlanks(self, rowNum):
        row = ""
        for btypeInt in self.sortedBondTypesAsInt:
            try:
                blength = self.bondLenghtsByBondType[btypeInt][rowNum]
            except:
                blength = "NA"
                # print(btypeInt, rowNum)
                # sys.exit()
            row += str(blength) + "\t"
        return row


def splitEntry(entry):
    tokens = entry.split(" ")
    return tokens[1], tokens[2]
