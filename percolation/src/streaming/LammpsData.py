#!/usr/bin/env python3

import sys
import re  # regex
import numpy as np

np.set_printoptions(edgeitems=20)

# Reading the system data file for system (lammps data file)
class LammpsData:
    nimpropers = 0  # treat impropers as dihedrals

    def __init__(self, dataFile, options):
        """ attributes """
        # instance variables: unique to each instance of Lammps
        # must be preceded by self.

        dataFileName = dataFile

        # num of atoms/bonds/angles/dihedrals
        self.numAtoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0

        # num of types
        self.natomTypes = 0
        self.nbondt = 0
        self.nanglet = 0
        self.ndihedralt = 0

        # box dimensions
        self.xlo = 0.0
        self.xhi = 0.0
        self.ylo = 0.0
        self.yhi = 0.0
        self.zlo = 0.0
        self.zhi = 0.0

        # box lengths
        self.boxLen = [0.0, 0.0, 0.0]

        # force field and atomic Masses

        # (for 0...natomTypes-1)
        self.atomTypeNames = []  # amber atom type names
        self.atomicMasses = []
        self.chemicalElements = [] # optional
        self.totalMass = 0
        self.paircoeffs = [[], []]  # sigma, epsilon

        # (for 0...nbondt-1)
        self.bondcoeffs = [[], []]  # force const, equil length

        # (for 0...nanglet-1)
        self.anglecoeffs = [[], []]  # force const, equil angle

        # (for 0...dihedralt-1)
        self.dihedralcoeffs = [[], [], []]  # K,d,n

        # per atom information (for 0...numAtoms-1) - numpy arrays
        self.charge = []  # partial atomic charges
        self.atomTypes = []  # lammps atom type numbers
        self.x = []  # coordinates

        # molecular topology - numpy arrays
        self.bonds = []  # (for 0...nbonds-1)
        self.angles = []  # (for 0...nangles-1)
        self.dihedrals = []  # (for 0...ndihedrals-1)

        # parsed flags - equal to 1 if section has been already parsed,
        # 0 otherwise
        # these sections are read only if needed

        self.parsedAtoms = 0
        self.parsedBonds = 0
        self.parsedAngles = 0
        self.parsedDihedrals = 0

        """ parse data file header + force field sections """

        with open(dataFile) as data:

            count = 0
            nsections = 16

            while count < nsections:

                l = data.readline()
                l = l.rstrip(" \n")  # remove trailing white spaces and EOL

                pos = l.find("atoms")
                if pos > -1:
                    self.numAtoms = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("atom t")
                if pos > -1:
                    self.natomTypes = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("bonds")
                if pos > -1:
                    self.nbonds = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("bond t")
                if pos > -1:
                    self.nbondt = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("angles")
                if pos > -1:
                    self.nangles = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("angle t")
                if pos > -1:
                    self.nanglet = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("dihedrals")
                if pos > -1:
                    self.ndihedrals = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("dihedral t")
                if pos > -1:
                    self.ndihedralt = int(l[0 : pos - 1])
                    count += 1
                    continue

                pos = l.find("xlo")
                if pos > -1:
                    l = l[0 : pos - 1]
                    l = l.split(" ")
                    self.xlo = float(l[0])
                    self.xhi = float(l[1])
                    count += 1
                    continue

                pos = l.find("ylo")
                if pos > -1:
                    l = l[0 : pos - 1]
                    l = l.split(" ")
                    self.ylo = float(l[0])
                    self.yhi = float(l[1])
                    count += 1
                    continue

                pos = l.find("zlo")
                if pos > -1:
                    l = l[0 : pos - 1]
                    l = l.split(" ")
                    self.zlo = float(l[0])
                    self.zhi = float(l[1])
                    count += 1
                    continue

                pos = l.find("Masses")
                if pos > -1:
                    self.readAtomicMassesAndChemicalElements(data)
                    count += 1
                    continue

                pos = l.find("Pair Coeffs")
                if pos > -1:
                    self.readPairCoeffs(data)
                    count += 1
                    continue

                pos = l.find("Bond Coeffs")
                if pos > -1:
                    self.readBondCoeffs(data)
                    count += 1
                    continue

                pos = l.find("Angle Coeffs")
                if pos > -1:
                    self.readAngleCoeffs(data)
                    count += 1
                    continue

                pos = l.find("Dihedral Coeffs")
                if pos > -1:
                    self.readDihedralCoeffs(data)
                    count += 1
                    continue

            if count < nsections:
                print(
                    "ERROR: one or more sections absent from header and/or",
                    "force field coeffs of data file. Quitting.",
                )
                sys.exit()

        self.boxLen[0] = self.xhi - self.xlo
        self.boxLen[1] = self.yhi - self.ylo
        self.boxLen[2] = self.zhi - self.zlo

        for op in options:

            if op == "atoms":
                self.readAtoms(dataFile)
                continue

            if op == "bonds":
                self.readBonds(dataFile)
                continue

            if op == "bondsLammps":
                self.readBondsLammpsStyle(dataFile)
                continue

            if op == "angles":
                self.readAngles(dataFile)
                continue

            if op == "dihedrals":
                self.readDihedrals(dataFile)
                continue

        print("\nLammps data file", dataFileName, "read in, with options:")
        for op in options:
            print(op)

    """ Parse methods """

    def readAtomicMassesAndChemicalElements(self, data):

        data.readline()

        for index in range(1, self.natomTypes + 1):
            l = data.readline()
            l = l.rstrip(" \n")  # remove trailing white spaces and EOL
            l = l.split(" ")
            self.atomicMasses.append(float(l[1]))
            try:
                self.chemicalElements.append(l[3])
            except:
                continue

    def readPairCoeffs(self, data):

        data.readline()

        for index in range(1, self.natomTypes + 1):
            l = data.readline()
            l = l.rstrip(" \n")
            l = l.split(" ")
            self.paircoeffs[0].append(float(l[1]))
            self.paircoeffs[1].append(float(l[2]))

    def readBondCoeffs(self, data):

        data.readline()

        for index in range(1, self.nbondt + 1):
            l = data.readline()
            l = l.rstrip(" \n")
            l = l.split(" ")
            self.bondcoeffs[0].append(float(l[1]))
            self.bondcoeffs[1].append(float(l[2]))

    def readAngleCoeffs(self, data):

        data.readline()

        for index in range(1, self.nanglet + 1):
            l = data.readline()
            l = l.rstrip(" \n")
            l = l.split(" ")
            self.anglecoeffs[0].append(float(l[1]))
            self.anglecoeffs[1].append(float(l[2]))

    def readDihedralCoeffs(self, data):

        data.readline()

        for index in range(1, self.ndihedralt + 1):
            l = data.readline()
            l = l.rstrip(" \n")
            l = l.split(" ")
            self.dihedralcoeffs[0].append(float(l[1]))
            self.dihedralcoeffs[1].append(int(l[2]))
            self.dihedralcoeffs[2].append(int(l[3]))

    def readAtoms(self, dataFile):

        self.atomTypes = np.zeros(self.numAtoms, dtype=int)
        self.charge = np.zeros(self.numAtoms)
        self.x = np.zeros((self.numAtoms, 3))

        with open(dataFile) as data:

            # since sections can come in any order, searching for the section
            # is the safest option (albeit not as efficient)

            l = data.readline()
            l = l.rstrip(" \n")
            while l.find("Atoms") == -1:
                l = data.readline()
                l = l.rstrip(" \n")
            data.readline()

            for index in range(self.numAtoms):
                l = data.readline()
                l = l.rstrip(" \n")
                line = l.split(" ")

                atomId = int(line[0])

                self.atomTypes[atomId - 1] = int(line[2])
                self.totalMass += self.atomicMasses[self.atomTypes[atomId - 1] - 1]
                self.charge[atomId - 1] = float(line[3])
                self.x[atomId - 1][0] = float(line[4])  # x
                self.x[atomId - 1][1] = float(line[5])  # y
                self.x[atomId - 1][2] = float(line[6])  # z

        self.parsedAtoms = 1

    def readBonds(self, dataFile):

        self.bonds = np.zeros((self.nbonds, 3), dtype=int)

        with open(dataFile) as data:

            l = data.readline()
            l = l.rstrip(" \n")
            while l.find("Bonds") == -1:
                l = data.readline()
                l = l.rstrip(" \n")
            data.readline()

            for index in range(self.nbonds):
                l = data.readline()
                l = l.rstrip(" \n")
                line = l.split(" ")

                bondId = int(line[0])

                self.bonds[bondId - 1][0] = int(line[1])  # bond type
                self.bonds[bondId - 1][1] = int(line[2])  # atom 1
                self.bonds[bondId - 1][2] = int(line[3])  # atom 2

        self.parsedBonds = 1

    def readBondsLammpsStyle(self, dataFile):

        self.bond_atom = {}  # key: atom id, value: list of atoms bonded to id
        self.bond_type = {}  # key: atom id, value: bond types of bonds involving id

        for i in range(self.nbonds):
            self.bond_atom[i] = []
            self.bond_type[i] = []

        with open(dataFile) as data:

            l = data.readline()
            l = l.rstrip(" \n")
            while l.find("Bonds") == -1:
                l = data.readline()
                l = l.rstrip(" \n")
            data.readline()

            for index in range(self.nbonds):
                l = data.readline()
                l = l.rstrip(" \n")
                line = l.split(" ")
                btype = int(line[1])
                id1 = int(line[2])
                id2 = int(line[3])
                self.bond_atom[id1 - 1].append(id2)
                self.bond_atom[id2 - 1].append(id1)
                self.bond_type[id1 - 1].append(btype)
                self.bond_type[id2 - 1].append(btype)

        self.parsedBondsLammps = 1

    def readAngles(self, dataFile):

        self.angles = np.zeros((self.nangles, 4), dtype=int)

        with open(dataFile) as data:

            l = data.readline()
            l = l.rstrip(" \n")
            while l.find("Angles") == -1:
                l = data.readline()
                l = l.rstrip(" \n")
            data.readline()

            for index in range(self.nangles):
                l = data.readline()
                l = l.rstrip(" \n")
                line = l.split(" ")

                angleId = int(line[0])

                self.angles[angleId - 1][0] = int(line[1])  # angle type
                self.angles[angleId - 1][1] = int(line[2])  # atom 1
                self.angles[angleId - 1][2] = int(line[3])  # atom 2
                self.angles[angleId - 1][3] = int(line[4])  # atom 3

        self.parsedAngles = 1

    def readDihedrals(self, dataFile):

        self.dihedrals = np.zeros((self.ndihedrals, 5), dtype=int)

        with open(dataFile) as data:

            l = data.readline()
            l = l.rstrip(" \n")
            while l.find("Dihedrals") == -1:
                l = data.readline()
                l = l.rstrip(" \n")
            data.readline()

            for index in range(self.ndihedrals):
                l = data.readline()
                l = l.rstrip(" \n")
                line = l.split(" ")

                dihedralId = int(line[0])

                self.dihedrals[dihedralId - 1][0] = int(line[1])  # dihedral type
                self.dihedrals[dihedralId - 1][1] = int(line[2])  # atom 1
                self.dihedrals[dihedralId - 1][2] = int(line[3])  # atom 2
                self.dihedrals[dihedralId - 1][3] = int(line[4])  # atom 3
                self.dihedrals[dihedralId - 1][4] = int(line[5])  # atom 4

        self.parsedDihedrals = 1

    """ Print methods """

    def printHeader(self):

        print(self.numAtoms, "atoms")
        print(self.nbonds, "bonds")
        print(self.nangles, "angles")
        print(self.ndihedrals, "dihedrals\n")

        print(self.natomTypes, "atom types")
        print(self.nbondt, "bond types")
        print(self.nanglet, "angle types")
        print(self.ndihedralt, "dihedral types\n")

        print(self.xlo, self.xhi, "xlo xhi")
        print(self.ylo, self.yhi, "ylo yhi")
        print(self.zlo, self.zhi, "zlo zhi")

    def printFF(self):

        print("\natomic Masses\n")
        for index in range(self.natomTypes):
            print(index + 1, self.atomicMasses[index])

        print("\nPair Coeffs\n")
        for index in range(self.natomTypes):
            print(index + 1, self.paircoeffs[0][index], self.paircoeffs[1][index])

        print("\nBond Coeffs\n")
        for index in range(self.nbondt):
            print(index + 1, self.bondcoeffs[0][index], self.bondcoeffs[1][index])

        print("\nAngle Coeffs\n")
        for index in range(self.nanglet):
            print(index + 1, self.anglecoeffs[0][index], self.anglecoeffs[1][index])

        print("\nDihedral Coeffs\n")
        for index in range(self.ndihedralt):
            print(
                index + 1,
                self.dihedralcoeffs[0][index],
                self.dihedralcoeffs[1][index],
                self.dihedralcoeffs[2][index],
            )

    def printAtoms(self):

        if not self.parsedAtoms:
            print("ERROR: no atoms in memory. Please call readAtoms()")
            return 1

        print("Atoms")
        for index in range(self.numAtoms):
            print(
                index + 1,
                "0",
                self.atomTypes[index],
                self.charge[index],
                self.x[index][0],
                self.x[index][1],
                self.x[index][2],
            )

    def printBonds(self):

        if not self.parsedBonds:
            print("ERROR: no bonds in memory. Please call readBonds()")
            return 1

        print("Bonds\n")
        for index in range(self.nbonds):
            print(
                index + 1,
                self.bonds[index][0],
                self.bonds[index][1],
                self.bonds[index][2],
            )

    def printAngles(self):

        if not self.parsedAngles:
            print("ERROR: no angles in memory. Please call readAngles()")
            return 1

        print("Angles\n")
        for index in range(self.nangles):
            print(
                index + 1,
                self.angles[index][0],
                self.angles[index][1],
                self.angles[index][2],
                self.angles[index][3],
            )

    def printDihedrals(self):

        if not self.parsedDihedrals:
            print("ERROR: no dihedrals in memory. Please call readDihedrals()")
            return 1

        print("Dihedrals\n")
        for index in range(self.ndihedrals):
            print(
                index + 1,
                self.dihedrals[index][0],
                self.dihedrals[index][1],
                self.dihedrals[index][2],
                self.dihedrals[index][3],
                self.dihedrals[index][4],
            )

    """ Analysis methods """

    # For each bond type, output mean, sd, min and max

    def bondStats(self):

        print("\nComputing bond length statistics\n")

        if not self.parsedAtoms:
            self.readAtoms(self.dataFileName)
        if not self.parsedBonds:
            self.readBonds(self.dataFileName)

        numBonds = np.zeros(self.nbondt, dtype=int)  # num of bonds of each
        # type, for 1...nbondt
        min = []
        for i in range(self.nbondt):
            min.append(200.0)
        min = np.array(min)  # min bond length of each type, for 1...nbondt

        shortestBondAtoms = np.zeros((self.nbondt, 2), dtype=int)

        max = np.zeros(self.nbondt)  # max bond length of each type,
        # for 1...nbondt

        longestBondAtoms = np.zeros((self.nbondt, 2), dtype=int)

        sd = np.zeros(self.nbondt)  # sd of bond lenghts of each type,
        # for 1...nbondt

        mean = np.zeros(self.nbondt)  # mean of bond lengths of each type,
        # for 1...nbondt

        sorted = np.zeros((self.nbondt, self.nbonds))  # multidim array to host
        # bonds sorted by type

        # sort bonds by type, compute min, max and sum
        for bond in self.bonds:

            btype = bond[0]
            dist = self.dist(bond[1], bond[2])

            mean[btype - 1] += dist
            sorted[btype - 1][numBonds[btype - 1]] = dist
            numBonds[btype - 1] += 1

            if dist < min[btype - 1]:
                min[btype - 1] = dist
                shortestBondAtoms[btype - 1][0] = bond[1]
                shortestBondAtoms[btype - 1][1] = bond[2]

            if dist > max[btype - 1]:
                max[btype - 1] = dist
                longestBondAtoms[btype - 1][0] = bond[1]
                longestBondAtoms[btype - 1][1] = bond[2]

        # compute mean and sd for each bond type
        for i in range(self.nbondt):

            if numBonds[i] == 0:
                print("No bonds of type", i + 1, "found.")
                continue

            mean[i] /= numBonds[i]

            for j in range(numBonds[i]):
                dist = sorted[i][j] - mean[i]
                sd[i] += dist * dist

            sd[i] /= numBonds[i] - 1
            sd[i] = np.sqrt(sd[i])

        # printing
        print("\n\t\t\tBond statistics:\n")
        print(
            "bond type  equilibrium  min  minIds  %(equilibrium-min)  mean",
            "  sd  max  maxIds  %(max-equilibrium)",
            sep="",
        )

        for i in range(self.nbondt):

            if numBonds[i] == 0:  # no bonds of this type found
                continue

            equilibrium = self.bondcoeffs[1][i]
            maxp = (max[i] - equilibrium) / equilibrium * 100
            minp = (equilibrium - min[i]) / equilibrium * 100

            print(
                i + 1,
                equilibrium,
                "{:0.4f}".format(min[i]),
                shortestBondAtoms[i],
                "{:0.1f}".format(minp),
                "{:0.4f}".format(mean[i]),
                "{:0.4f}".format(sd[i]),
                "{:0.4f}".format(max[i]),
                longestBondAtoms[i],
                "{:0.1f}".format(maxp),
                sep="\t",
            )

    # compute the Euclidean distance between two atoms
    def dist(self, atom1, atom2):

        coord1 = [self.x[atom1 - 1][0], self.x[atom1 - 1][1], self.x[atom1 - 1][2]]
        coord2 = [self.x[atom2 - 1][0], self.x[atom2 - 1][1], self.x[atom2 - 1][2]]
        l = [0.0, 0.0, 0.0]

        for i in [0, 1, 2]:

            # swap so that x2 > x1
            if coord1[i] > coord2[i]:
                tmp = coord1[i]
                coord1[i] = coord2[i]
                coord2[i] = tmp

            l[i] = coord2[i] - coord1[i]
            # check if bond straddles a periodic boundary
            # heuristic: check if bond is way too long
            if l[i] > self.boxLen[i] / 2:
                l[i] = coord1[i] - coord2[i] + self.boxLen[i]

        return np.sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2])
