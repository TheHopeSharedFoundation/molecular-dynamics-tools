#/usr/local/bin/pymol
#import pymol
#from pymol import cmd
#from pymol import stored
import numpy as np
import sklearn
from sklearn.decomposition import PCA
from scipy.sparse import coo_matrix
import scoria
import mdtraj as md
from sklearn.decomposition import TruncatedSVD
from scipy import sparse as sp
from scipy.sparse import csr_matrix
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
import argparse






traj = md.load('/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.xtc', top='/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.pdb')

####### Load the aggregate pocket and Ne-atom platform

pocket = scoria.Molecule()
platform = scoria.Molecule()
pocket.load_pdb_into("/Users/benjaminsamudio/Desktop/EGFR_inactiveState_pocketOfInterest.pdb") # This is the "aggregate" pocket which comes from the combination of all cosolvent trajectories
platform.load_pdb_into("/Users/benjaminsamudio/Desktop/Ne_atom_grid_0p50_angstrom_spacing_2d_testing.pdb") # This is Ne atom 2D platform that moves across the area of the aggregate pocket 

test_array = []

#single_frame = md.load_frame('/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.xtc',index=args.frame,top='/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.pdb')

translation_vector = np.array([0,0,-0.5]) # This is the vector that describes translation across the area of the aggregate pocket
reset_vector = np.array([0,0,25]) # This is the vector that describes the reseting of the position of the platform
#for frame in range(1,int(traj.n_frames)+1): # Each frame takes about seconds to process.  In a typical trajectory, there are about 42000 frames.  
#for frame in range(1,10):
#single_frame = traj[frame]
	
#single_frame.save_pdb('single_frame.pdb') # append a unique identifier 

with open("grid_based_fingerprints-TEST.csv","w") as fingerprint_out:
	fingerprint_out.write(f"frame,size,fingerprint\n")
	for frame in range(0,int(traj.n_frames)): #<--------------------- changed int(traj.n_frames) to a test integer value
		single_frame = traj[frame]
		single_frame.save('single_frame.pdb')
		protein = scoria.Molecule()
		protein.load_pdb_into("single_frame.pdb")
		print(f"Operating on frame: {frame}")
		row = []
		column = []
		value = []
		output_string = ""
		fingerprint_size = 0
		out_mol = scoria.Molecule()
		for i in range(0,50):
			platform.translate_molecule(translation_vector)
			platform_clash_protein = platform.select_close_atoms_from_different_molecules(protein,1.8)
			platform_clash_pocket = platform.select_close_atoms_from_different_molecules(pocket,3)
			platform_defined_pocket = list(set(platform_clash_pocket[0]) - set(platform_clash_protein[0]))
			if len(platform_defined_pocket) > 0:
				print(platform_defined_pocket)
				#test_mol = scoria.Molecule()
				#test_mol = platform.get_molecule_from_selection(platform_defined_pocket)
				#test_file_name = "TestFile_" + str(frame) + "_" + str(i) + ".pdb"
				#test_mol.save_pdb(test_file_name)
				for j in range(0,len(platform_defined_pocket)):
					fingerprint_size += 1
					#output_string = output_string + str(i)  + "#" + str(platform_defined_pocket[j]) + " "
					output_string = output_string  + str(platform_defined_pocket[j]) + " "
		fingerprint_out.write(f"{frame},{fingerprint_size},{output_string}\n")
		platform.translate_molecule(reset_vector)

