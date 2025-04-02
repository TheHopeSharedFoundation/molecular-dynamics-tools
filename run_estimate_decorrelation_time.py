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
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import KMeans
import argparse
import pandas as pd
import pandas as pd
import plotly.express as px
import random
from collections import Counter

####################################  Testing start

# The following is based on: https://pmc.ncbi.nlm.nih.gov/articles/PMC2538559/

traj = md.load("/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.dcd",top="/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.pdb")

print(f"The number of frames is: {traj.n_frames}")

df = pd.DataFrame(np.zeros((int(traj.n_frames), 2)), columns=['structure_bin_index', 'rmsd'])
df['structure_bin_index'] = -1
df['rmsd'] = -1

print(df.head())

number_references = 10
subsample_size = 2
picoseconds_per_frame = 10
fractional_occupancy = 1 / number_references
number_frames = int(traj.n_frames)
block_size = int(number_frames * fractional_occupancy)
iid = ((1/number_references-1/number_references**2)/subsample_size) * ((number_frames-subsample_size)/(number_frames-1)) 

print(f"This is iid: {iid} {number_frames}")

############################### Initialize the bin array

protein_indices = []
protein_indices = traj.topology.select("protein")

for bin_index in range(0,number_references):
	# Filter DataFrame
	filtered_df = df[df['structure_bin_index'] == -1]

	random_index = random.randint(0, len(filtered_df))
	print(f"Random row: {random_index} bin_index: {bin_index}\n")
	rmsds = md.rmsd(traj, traj, random_index, atom_indices=protein_indices) 
	for index in range(0,len(filtered_df)):
		#similarity = jaccard_similarity(filtered_df['fingerprint'].iloc[random_index],filtered_df['fingerprint'].iloc[index])
		similarity = rmsds[index]
		filtered_df['rmsd'].iloc[index] = similarity
		filtered_df['structure_bin_index'].iloc[index] = bin_index

	filtered_df.sort_values(by='rmsd',inplace=True)
	filtered_df = filtered_df.head(block_size)
	print(filtered_df)
	df.update(filtered_df)
	
print(df)

###############################

observed_list = []
frame_list = []
separations_per_trajectory = 0
number_subsamples = 0

for frame_separation in range(0,number_frames): #<----------------- Loop through each frame of the trajectory.  This is the same as the frame separation.
	separations_per_trajectory = int(number_frames / (frame_separation+1))
	number_subsamples = int(separations_per_trajectory / subsample_size)
	if number_subsamples >= 2:
		bin_frame = []
		bin_instance_list = [] #<------------------ This list simply captures the index of the bin
		subsample_fill_count = 0 #<------- This is used to test that the subsample size has reached that of the set/specified subsample size
		whole_trajectory_index = 0
		for subsample_count in range(0,number_subsamples):
			bin_instance_list.append(int(df['structure_bin_index'].iloc[whole_trajectory_index]))
			whole_trajectory_index = whole_trajectory_index + (frame_separation + 1)
			subsample_fill_count += 1
			if subsample_fill_count == subsample_size: # <------------ Test whether the size of the subsample matches the set/specified subsample size
				bin_subsample = [0] * number_references #<----------------- Create a temporary list of zeros in which the element index = the bin index.
				bin_frequencies = Counter(bin_instance_list)
				for key, value in bin_frequencies.items():
					bin_subsample[key] = value 
				bin_frame.append(bin_subsample)
				#print(bin_subsample)
				bin_instance_list = []
				subsample_fill_count = 0
		#print(f"frame: {frame_separation} length {len(bin_frame)}")
		#print(bin_frame)
		real_variance = 0
		observed_variance = 0
		if len(bin_frame) >= 1:
			real_variance_avg_sum = 0
			real_variance_avg = 0
			for bin_index in range(0,number_references):
				real_variance_sum = 0
				real_variance  = 0
				for subsample_index in range(0,len(bin_frame)):
					#print(f"{bin_index},{bin_frame[subsample_index][bin_index]}")
					real_variance_sum = real_variance_sum + (bin_frame[subsample_index][bin_index]/subsample_size - 1/number_references)**2
				real_variance = real_variance_sum / len(bin_frame)
				#print(f"{real_variance}")
				real_variance_avg_sum = real_variance_avg_sum + real_variance 
			real_variance_avg = real_variance_avg_sum / number_references
			#print(f"{real_variance_avg_sum},{real_variance_avg}")
			observed_variance = real_variance_avg / iid
			print(observed_variance)

