############################################################################
# Copyright 2024 The Hope Shared Foundation                                #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
# http://www.apache.org/licenses/LICENSE-2.0                               #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
#                                                                          #
# The Hope Shared Foundation, September 2024                               #
# https://hopesharedfoundation.org/                                        #
# Towards Alleviating Suffering                                            #
############################################################################

import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.analysis import rms, align
import pandas as pd
import numpy as np
import os
import argparse
from datetime import datetime
from pathlib import Path

vdwradii = {'Aa': 1.85, 'Ag': 1.72, 'Al': 1.84, 'Ar': 1.88, 'At': 2.02, 'Au': 1.66, 'B': 1.92, 'Ba': 2.68, 'Be': 1.53, 'Bi': 2.07, 'Br': 1.85, 'C': 1.7, 'Ca': 2.31, 'Cd': 1.58, 'Cl': 1.75, 'Cs': 3.43, 'Cu': 1.4, 'F': 1.47, 'Fr': 3.48, 'Ga': 1.87, 'Ge': 2.11, 'H': 1.1, 'He': 1.4, 'Hh': 1.55, 'I': 1.98, 'In': 1.93, 'K': 2.75, 'Kr': 2.02, 'Li': 1.82, 'Mg': 1.73, 'N': 1.55, 'Na': 2.27, 'Ne': 1.54, 'Ni': 1.63, 'O': 1.52, 'P': 1.8, 'Pb': 2.02, 'Pd': 1.63, 'Po': 1.97, 'Pt': 1.75, 'Ra': 2.83, 'Rn': 2.2, 'Rr': 3.03, 'S': 1.8, 'Sb': 2.06, 'Se': 1.9, 'Si': 2.1, 'Sn': 2.17, 'Sr': 2.49, 'Te': 2.06, 'Tl': 1.96, 'U': 1.86, 'Xe': 2.16, 'Zn': 1.39}

# Source: https://docs.mdanalysis.org/2.0.0/documentation_pages/topology/tables.html#MDAnalysis.topology.tables.vdwradii

####################### References

# https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/

####################### Set a path to the trajectory file (DCD format).  A PDB structure is input to OpenMM.  To this structure, OpenMM adds water molecules.  This structure is considered the "initial" structure.  Set a path to the initial structure (PDB format).

parser = argparse.ArgumentParser(prog='Transform post-MD trajectory coordinates',description='This program transforms coordinates from an MD trajectory',epilog='Copyright 2024 The Hope Shared Foundation')

parser.add_argument('--trajectory',help="This is a molecular dynamics trajectory file.  This argument is a required argument.",required=True)
parser.add_argument('--topology',help="This is a topology structure in the protein data bank (PDB) file format.  This is a required argument.",required=True)
parser.add_argument('--reference',help="This is a reference structure in the protein data bank (PDB) file format to which trajectory structures will be translated and rotated.  This is a required argument.",required=True)
parser.add_argument('--skip',default="NULL",help="This is the trajectory frame skip size.  The default is NO skip.",required=False)
parser.add_argument('--alpha',default="NULL",help="This is a flag to calculate alpha-spheres and ligand coordinates or not.  The default is NO calculation.  To enable calculation, type yes after this argument.",required=False)
args = parser.parse_args()

TRJ = args.trajectory
TOP = args.topology
REF = args.reference

######################## Read the trajectory and initial structure into memory

if args.skip == "NULL":
    u = mda.Universe(TOP,TRJ,guess_bonds=True,vdwradii=vdwradii)
    ref = mda.Universe(REF,guess_bonds=True,vdwradii=vdwradii)
    u.transfer_to_memory(verbose=True)
    ref.transfer_to_memory(verbose=True)
else:
    skip_size = int(args.skip)
    u = mda.Universe(TOP,TRJ,guess_bonds=True,vdwradii=vdwradii)
    ref = mda.Universe(REF,guess_bonds=True,vdwradii=vdwradii)
    u.transfer_to_memory(verbose=True,step=skip_size)
    ref.transfer_to_memory(verbose=True,step=skip_size)

print(f"Number of atoms: {u}")
print(f"Number of frames: {len(u.trajectory)}")

####################### Define sets of atoms.  Benzene molecules have a residue name of "UNK".

reference = ref.select_atoms("protein") + ref.select_atoms("resname UNK")
reference_CA = ref.select_atoms("protein and name CA")
protein_CA = u.select_atoms("protein and name CA")
solute = u.select_atoms("protein") + u.select_atoms("resname UNK")
protein = u.select_atoms("protein")
ligands = u.select_atoms("resname UNK")


####################### Apply structural transformations: 1) keep complexes together with "NoJump", 2) center molecules in a box, 3) wrap molecules into the box, 4) fit, rotate, and translate the protein to a reference structure

# originally after "noJump": transformations.unwrap(solute,parallelizable=True)
# originally before "unwrap": transformations.nojump.NoJump(parallelizable=True)


#workflow = (transformations.unwrap(solute,parallelizable=True),transformations.center_in_box(protein,center='mass'),transformations.wrap(solute,compound='fragments',parallelizable=True),transformations.fit_rot_trans(protein,reference,parallelizable=True))
#u.trajectory.add_transformations(*workflow)



workflow = [transformations.unwrap(solute,parallelizable=True),transformations.center_in_box(protein, center='mass',parallelizable=True),transformations.wrap(solute, compound='fragments',parallelizable=True),transformations.fit_rot_trans(protein_CA,reference_CA,parallelizable=True)]
u.trajectory.add_transformations(*workflow)

####################### Write the structurally transformed trajectory to a new DCD file


output_trajectory_dcd = args.trajectory + "_transformed_skip_" + str(args.skip) + ".dcd"
output_trajectory_xtc = args.trajectory + "_transformed_skip_" + str(args.skip) + ".xtc"
output_topology = args.trajectory + "_transformed_skip_" + str(args.skip) + ".pdb"
solute.write(output_trajectory_dcd,frames=u.trajectory[:])
solute.write(output_trajectory_xtc,frames=u.trajectory[:])
solute.write(output_topology,frames=u.trajectory[1:2])

if args.alpha == "yes":

    path = Path(args.trajectory)
    file_stem = path.stem
    atoms_O_file = file_stem + "_transformed_skip_" + str(args.skip) + "_oxygenAtoms.txt"
    atoms_N_file = file_stem + "_transformed_skip_" + str(args.skip) + "_nitrogenAtoms.txt"
    atoms_C_file = file_stem + "_transformed_skip_" + str(args.skip) + "_carbonAtoms.txt"


    with open(atoms_O_file,'w') as atoms_O_output:
        with open(atoms_N_file,'w') as atoms_N_output:
            with open(atoms_C_file,'w') as atoms_C_output:

                atoms_O = u.select_atoms("resname UNK and element O")
                atoms_N = u.select_atoms("resname UNK and element N")
                atoms_C = u.select_atoms("resname UNK and element C")
                
                output_index = 0
                dir_name_pqr = file_stem + "_pqr_files"
                cmd_rm = "rm -r " + dir_name_pqr
                os.system(cmd_rm)
                cmd_mkdir = "mkdir " + dir_name_pqr
                os.system(cmd_mkdir)
                dir_name_pdb = file_stem + "_pdb_files"
                cmd_rm = "rm -r " + dir_name_pdb
                os.system(cmd_rm)
                cmd_mkdir = "mkdir " + dir_name_pdb
                os.system(cmd_mkdir)
                
                # Get coordinates for atoms from trajectories using MDanalysis: https://stackoverflow.com/questions/66660281/using-mdanalysis-to-extract-coordinates-in-an-array-from-pdb

                for ts in u.trajectory:
                    atoms_O_coordinates = atoms_O.atoms.positions
                    atoms_N_coordinates = atoms_N.atoms.positions
                    atoms_C_coordinates = atoms_C.atoms.positions
                    print(atoms_O_coordinates, file=atoms_O_output)
                    print(atoms_N_coordinates, file=atoms_N_output)
                    print(atoms_C_coordinates, file=atoms_C_output)
                    filename = file_stem + '_' + str(output_index) + '.pdb'
                    filebase = file_stem + '_' + str(output_index)
                    solute.write(filename,frames=u.trajectory[output_index:output_index+1])
                    cmd = "fpocket -f " + filename
                    os.system(cmd)
                    cmd = "mv ./" + filebase + "_out/" + filebase + "_pockets.pqr" + " ./" + dir_name_pqr
                    os.system(cmd)
                    cmd = "mv ./" + filename + " ./" + dir_name_pdb
                    os.system(cmd)
                    cmd = "rm -r " + filebase + "_out/"
                    os.system(cmd)
                    output_index += 1
