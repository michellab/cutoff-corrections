#
# Computes correction terms to LJ energies
#
import sys,os
from Sire.IO import *
from Sire.Mol import *
from Sire.Units import *
from Sire.MM import *


cutoff = 10.0 # Angstrom 
# List of (sig,eps) params for solvent sites
# How do we handle multiple VDW sites?
solvent_LJ_params = [ (1.0, 0.0) ]#Angstrom, kcal/mol
rho_solvent_model = 1.0# bulk density of model solvent (UNITS?)

solvent_residues = ["WAT","ZBK"]

def SplitSoluteSolvent(system):
    molecules = system.molecules()
    mol_numbers = molecules.molNums()
    solutes = MoleculeGroup("solutes")
    solvent = MoleculeGroup("solvent")
    for molnum in mol_numbers:
        mol = molecules.molecule(molnum).molecule()
        res0 = mol.residues()[0]
        if res0.name().value() in solvent_residues:
            solvent.add(mol)
        else:
            solutes.add(mol)
    return solutes, solvent

if __name__ == '__main__':
    # Step 1  Load a crd/top
    amber = Amber()
    system, space = amber.readCrdTop(crd_file, top_file)
    # Now filter out solvent molecules
    solutes, solvent = SplitSoluteSolvent(system)
