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

for frame in range(1,40000,1000):
	single_frame = traj[frame]
	single_frame.save('single_frame.pdb')
	protein = scoria.Molecule()
	protein.load_pdb_into("single_frame.pdb")
	print(f"Operating on frame: {frame}")
	row = []
	column = []
	value = []
	output_string = ""
	for i in range(0,50):
		platform.translate_molecule(translation_vector)
		platform_clash_protein = platform.select_close_atoms_from_different_molecules(protein,1.8)
		platform_avoid_protein = platform.invert_selection(platform_clash_protein[0])
		new_mol = scoria.Molecule()
		new_mol =  platform.get_molecule_from_selection(platform_avoid_protein)
		close_reference_pocket = new_mol.select_close_atoms_from_different_molecules(pocket,0.6)
		if len(close_reference_pocket[0]) > 0:
			for j in list(close_reference_pocket[0]):
				row.append(int(frame))
				column.append(str(i*1600+j))
				output_string = output_string + str(i*1600+j) + " "
				value.append(1)
	print(row)
	print(column)
	print(value)
	test_array.append(output_string)
	platform.translate_molecule(reset_vector)



#print(test_array)

#tfidf_vectorizer = TfidfVectorizer(max_features=5000)
#X_tfidf = tfidf_vectorizer.fit_transform(test_array)

#print(X_tfidf)

text_features = HashingVectorizer(analyzer="word").fit_transform(test_array)
#print(text_features)
clusterer = KMeans(n_clusters=2)
clusters = clusterer.fit_predict(text_features)
X_pca = PCA(n_components=2).fit_transform(text_features.toarray())

pca_x_list = []
pca_y_list = []

for i in range(0,len(X_pca)):
	#print(f"{X_pca[i][0]},{X_pca[i][1]},{clusters[i]}")
	pca_x_list.append(X_pca[i][0])
	pca_y_list.append(X_pca[i][1])


##############################

xedges = []
yedges = []


for increment in range(-100,100):
        xedges.append(increment*0.01)
        yedges.append(increment*0.01)

print(pca_x_list)
print(pca_y_list)

H = np.histogram2d(pca_x_list,pca_y_list, bins=(xedges,yedges))

#print(H)

H_max = np.max(H[0])

#print(H[0],H_max)

for i in range(0,len(H[0])):
	for j in range(0,len(H[0])):
		print(H[0][i][j])
		if H[0][i][j] != 0:
			pmf = -8.314 * 1/4.184 * 1/1000  * 300 * np.log(int(H[0][i][j])/int(H_max))
			print(pmf)

