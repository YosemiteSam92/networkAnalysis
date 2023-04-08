from lammpsData import LammpsData
from lammpsDataPrinter import LammpsDataPrinter as lprinter
from os import chdir, path

# change working directory to script location
chdir(path.dirname(path.realpath(__file__)))

dataFiles = [
    "data.equilibriate_298K_large1_127ns",
    "data.equilibriate_298K_large2_152ns",
    "data.equilibriate_298K_large3_152ns",
    "data.equilibriate_298K_large4_152ns",
    "data.equilibriate_298K_large5_107ns"
]

for dataFile in dataFiles:
    data = LammpsData("../lammpsData_annealedEquilibrated/"+dataFile)
    data.removeBondsAcrossPbcs()
    data.removeBondsByMass("1.008")
    dataPrinter = lprinter(data)
    dataPrinter.printDataFile("../lammpsData_networkAnalysis/"
                                +dataFile+"_noPbcBonds_noBondsWithH")