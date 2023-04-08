#!/usr/bin/env python3

import sys, resource
import argparse

if __name__ == "__main__":
    import sys
    from os.path import dirname, abspath

    sys.path.append(abspath(dirname(__file__)))

    from mols.molecularGraphFromMol2 import MolecularGraphFromMol2
    from mols.mol2Molecule import Mol2Molecule

    parser = argparse.ArgumentParser(
        description="Script to test percolation detection on sample molecules"
    )
    parser.add_argument(
        "molecule_index",
        type=int,
        default=15,
        help="Index of sample molecule to check.",
    )

    args = parser.parse_args()

    moleculeName = "molecule" + str(args.molecule_index)
    fileName = "molecule_data/" + moleculeName + ".mol2"
    mol2Molecule = Mol2Molecule(fileName)
    molGraph = MolecularGraphFromMol2(mol2Molecule)

    molGraph.isAnyMoleculePercolating()
    
    print("")
    print("Tested molecule:", moleculeName)
    print("Atoms in pbc components:", molGraph.atomsInPbcComponent)
    print("Number of percolating molecules:", molGraph.numPercolatingMolecules)
    
    molGraph.graph.dump(moleculeName + ".dot", True)
