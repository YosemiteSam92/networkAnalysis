#!/usr/bin/env python3

import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__)))
from streaming.LammpsData import LammpsData

dataFile = "lammpsData/data.Min"

lammps = LammpsData(dataFile, ["atoms", "bondsLammps"])
# lammps.countNumOfMolecules()
lammps.computeMolecularMasses()

# options = ["atoms","angles"]
# lammps2 = Lammps(dataFile,options)

# lammps.bondStats()
