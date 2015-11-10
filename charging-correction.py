#
# Computes correction terms to charging free energies
#
#
# DF_corr = DF_pol + DF_psum + DF_exc + DF_dir
#
import sys,os
from Sire.IO import *
from Sire.Mol import *
from Sire.Units import *
from Sire.MM import *
from Sire.Maths import *

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

def PoissonPBC(binary, solutes, space, cutoff, dielectric):
    
    space_x = space.dimensions()[0]/10.0 # nm
    space_y = space.dimensions()[1]/10.0 #
    space_z = space.dimensions()[2]/10.0 #
    nions = 0
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        nions += mol.nAtoms()
    dielec = dielectric
    cut = cutoff/10.0

    infilepart1 = """GRID
100 100 100
%8.5f %8.5f %8.5f
4
END
ITERATION
200 1.5 0.3 -0.001
0 50 50 50 1
END
ELECTRO
%s %8.5f
""" % (space_x, space_y, space_z, nions, dielec)

    infilepart2 = ""
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        atoms = mol.atoms()
        for atom in atoms:
            charge = atom.property("charge").value()
            sigma = atom.property("LJ").sigma().value()
            radius = 0.5*(sigma*2**(1/6.))/10.0 # to nm. Is this a decently good approximation?
            coords = atom.property("coordinates")/10.0 # to nm 
            line = "%8.5f %8.5f %8.5f %8.5f %8.5f\n" % (charge, radius, coords[0], coords[1], coords[2])
            #import pdb; pdb.set_trace()
            infilepart2 += line
    infilepart2 += "END\n" 

    infilepart3 = """BOUNDARY
3
%s
END
""" % (cut)
    infile = infilepart1 + infilepart2 + infilepart3
    #lines = infile.split('\n')
    #for line in lines:
    #    print (line)
    #import pdb; pdb.set_trace()

    wstream = open('infile','w')
    wstream.write(infile)
    wstream.close()

    cmd = "%s > PB.out" % binary
    #os.system(cmd)

    cmd = "tail -1 PB.out > temp"
    os.system(cmd)
    rstream = open('temp','r')
    buffer = rstream.readlines()
    elems = buffer[0].split()
    nrg = float(elems[1])

    DF_BA_PBC = nrg * kJ_per_mol

    return DF_BA_PBC

def PoissonNP(binary, solutes, dielectric):
    # Here write input file for apbs...
    DF_CB_NP = 0.0 * kcal_per_mol
    nions = 0
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()


    xmlinpart1 = """# Only the x, y, z, charge, and radii fields are read.
<ion-example>
  <residue>
     <resName>ION</resName>
     <chainID>A</chainID>
     <resSeq>1</resSeq>
"""

    xmlinpart2 = ""

    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        nions += mol.nAtoms()
        atoms = mol.atoms()
        for atom in atoms:
            name = atom.name().value()
            charge = atom.property("charge").value()
            sigma = atom.property("LJ").sigma().value()
            radius = 0.5*(sigma*2**(1/6.)) # Ang. Is this a decently good approximation?
            coords = atom.property("coordinates") # APBS uses angstroms
            xmlinpart2 += "     <atom>\n"
            xmlinpart2 += "        <serial>ATOM</serial>\n"
            xmlinpart2 += "        <name>%s</name>\n" % name
            xmlinpart2 += "        <x>%s</x>\n" % coords[0]
            xmlinpart2 += "        <y>%s</y>\n" % coords[1]
            xmlinpart2 += "        <z>%s</z>\n" % coords[2]
            xmlinpart2 += "        <charge>%s</charge>\n" % charge
            xmlinpart2 += "        <radius>%s</radius>\n" % radius
            xmlinpart2 += "     </atom>\n"

    xmlinpart3 = """  </residue>
</ion-example>
"""

    xmlin = xmlinpart1 + xmlinpart2 + xmlinpart3

    wstream = open("apbssystem.xml","w")
    wstream.write(xmlin)
    wstream.close()

    apbsin="""#############################################################################
### BORN ION SOLVATION ENERGY
### $Id$
###
### Please see APBS documentation (http://apbs.sourceforge.net/doc/) for 
### input file sytax.
#############################################################################

# READ IN MOLECULES
read                                               
    mol xml apbssystem.xml
end

# COMPUTE POTENTIAL FOR SOLVATED STATE
elec name solvated
    mg-auto     
    dime 65 65 65
    cglen 50 50 50
    fglen 12 12 12
    fgcent mol 1
    cgcent mol 1
    mol 1
    lpbe
    bcfl mdh
    pdie 1.0
    sdie %s
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy total
    calcforce no
    write pot gz potential
    # write pot dx potential
    # write charge dx charge
end

# COMPUTE POTENTIAL FOR REFERENCE STATE
elec name reference
    mg-auto
    dime 65 65 65
    cglen 50 50 50
    fglen 12 12 12
    fgcent mol 1
    cgcent mol 1
    mol 1
    lpbe
    bcfl mdh
    pdie 1.0
    sdie 1.0
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy total
    calcforce no
end

# COMBINE TO GIVE SOLVATION ENERGY
print elecEnergy solvated - reference end

quit
""" % dielectric

    wstream = open("apbs.in","w")
    wstream.write(apbsin)
    wstream.close()

    cmd = "%s apbs.in 1> apbs.out 2> apbs.err" % (binary)
    #os.system(cmd)

    rstream = open("apbs.out","r")
    buffer = rstream.readlines()

    nrg = 0.0
    nrgfound = False
    for line in buffer:
        if line.find("Global net ELEC" ) > 0:
            elems = line.split()
            #print (line)
            #print (elems)
            nrg = float(elems[5])
            nrgfound = True
            break
    if not nrgfound:
        print ("ERROR. Could not find electrostatic energy in apbs.out file. ABORT")
        sys.exit(-1)
    DF_CB_NP = nrg * kJ_per_mol

    return DF_CB_NP

def SummationCorrection(solutes, solvent, space, rho_solvent_model,
                        quadrupole_trace, eps_solvent_model, BAcutoff):
    nrg = 0.0
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()

    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        #import pdb; pdb.set_trace()
        atoms = mol.atoms()
        mol_charge = 0.0
        mol_com = [0.0, 0.0, 0.0]
        mol_mass = 0.0
        for atom in atoms:
            name = atom.name().value()
            charge = atom.property("charge").value()
            mass = atom.property("mass").value()
            coords = atom.property("coordinates")
            # What about PBC ? 
            mol_com[0] += coords[0]*mass
            mol_com[1] += coords[1]*mass
            mol_com[2] += coords[2]*mass
            mol_charge += charge
            mol_mass += mass
        for x in range(0,3):
            mol_com[x] = mol_com[x]/mol_mass
        #print (mol_com)
        #print (mol_charge)
        # JM 11/15 This code only deals with the first solute molecule at present
        break

    nsolv = 0
    solv_mols = solvent.molecules()
    solvmolnums = solv_mols.molNums()
    for molnum in solvmolnums:
        solvent = solv_mols.molecule(molnum).molecule()
        # Find com
        solv_com = [0.0, 0.0, 0.0]
        solv_atoms = solvent.atoms()
        solv_mass = 0.0
        for solv_atom in solv_atoms:
            coords = solv_atom.property("coordinates")
            mass = solv_atom.property("mass").value()
            solv_mass += mass
            for i in range(0,3):
                solv_com[i] += coords[i]*mass
        for i in range(0,3):
            solv_com[i] /= solv_mass
        #print (mol_com)
        #print (solv_com)
        d = space.calcDist(Vector(mol_com), Vector(solv_com))
        #print (d)
        if d < BAcutoff:
            nsolv += 1
    print (nsolv)
    ONE_OVER_6PI_EPS0 = 290.98622868361923
    nrg = -ONE_OVER_6PI_EPS0 * mol_charge * quadrupole_trace *\
        ( ( (2*(eps_solvent_model-1) / (2*eps_solvent_model+1) )*\
              nsolv/(4*pi*(BAcutoff/10.0)**3/3.0) ) +\
              (3/(2*eps_solvent_model)))
    # !!! MIissing (3/(2*eps_solvent_model+1)) phi_ODL term
    #sys.exit(-1)
    DF_PSUM = nrg * kJ_per_mol
    return DF_PSUM

######## THIS WILL EVENTUALLY BE PASSED AS INPUT FROM THE COMMAND LINE ##
top_file = "30x30x30-water/solvated.parm7"
crd_file = "30x30x30-water/md00004.rst7"

PoissonPBCSolverBin = "/home/julien/software/thirdparty/FFT_PB_SOLVER_FOR_CUTOFF+CORR/pb_generalT"
PoissonNPSolverBin = "/home/julien/software/thirdparty/APBS-1.4.1-binary/bin/apbs"

BAcutoff = 10.0# Value of the atom based Barker-Watt reaction field cutoff in Angstroms
eps_solvent_model = 78.4# dielectric constant of model solvent in MD simulations with PBC
eps_solvent_np = 78.4 # experimental dielectric constant of solvent in NP conditions
rho_solvent_model = 1.0# bulk density of model solvent for LJ corrections
quadrupole_trace = 0.0082# e.nm^2 !!!SPC needs adjustment
solvent_residues = ["WAT","ZBK"]
##########################################################################
if __name__ == '__main__':
    # Step 1  Load a crd/top
    amber = Amber()
    system, space = amber.readCrdTop(crd_file, top_file)
    # Now filter out solvent molecules
    solutes, solvent = SplitSoluteSolvent(system)
    # Step 2  Compute DF_pol
    # step 2b Use PH'solver to compute DF^{BA}_{PBC}
    DF_BA_PBC = PoissonPBC(PoissonPBCSolverBin ,solutes, space, BAcutoff, eps_solvent_model)
    print (DF_BA_PBC)
    # step 2a Use APBS to compute DF^{CB}_{infinite}
    DF_CB_NP = PoissonNP(PoissonNPSolverBin, solutes, eps_solvent_np)
    print (DF_CB_NP)
    #sys.exit(-1)
    DF_POL = DF_CB_NP - DF_BA_PBC
    print ("DF_POL IS %s " % DF_POL)
    # step 3  Compute psum
    DF_PSUM = SummationCorrection(solutes, solvent, space, rho_solvent_model,\
                                      quadrupole_trace, eps_solvent_model, \
                                      BAcutoff)
    print ("DF_PSUM IS %s" % DF_PSUM)
    sys.exit(-1)
    # step 3a Find number of solvent particles within Rcutoff distance
    # of the centre of mass of the discharged molecule
    # step 3b evaluate psum based on Kastenholz formulas
    # step 4 Compute DF_exc (Q? Do we need this ? )
    # how do we handle intramolecular electrostatics non-bonded in solutes?
    # might be we do not need it if used BA in vacuum leg (?)
    # only applies to polyatomic ions.
    # step 5 Compute DF_dir
    # must understand how to implement this term
