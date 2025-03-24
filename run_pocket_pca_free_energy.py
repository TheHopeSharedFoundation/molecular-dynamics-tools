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
import pandas as pd
import pandas as pd
import plotly.express as px
import random
from collections import Counter


def jaccard_similarity(doc1, doc2):
	words_doc1 = set(doc1.split())
	words_doc2 = set(doc2.split())
	intersection = words_doc1.intersection(words_doc2)
	union = words_doc1.union(words_doc2)
	return len(intersection) / len(union)

df = pd.read_csv("/Users/benjaminsamudio/grid_based_fingerprints.csv")

####################################  Testing start

# The following is based on: https://pmc.ncbi.nlm.nih.gov/articles/PMC2538559/

number_references = 10
number_subsamples = 10
picoseconds_per_frame = 10
fractional_occupancy = 1 / number_references
number_frames = len(df)
block_size = int(number_frames * fractional_occupancy)

#df =df.head(block_size)


# Initialize the bin array

df['jaccard_similarity'] = 0
df['structure_bin_index'] = -1

for bin_index in range(0,number_references):
	# Filter DataFrame
	filtered_df = df[df['structure_bin_index'] == -1]

	random_index = random.randint(0, len(filtered_df))
	print(f"Random row: {random_index} bin_index: {bin_index}\n")

	for index in range(0,len(filtered_df)):
		similarity = jaccard_similarity(filtered_df['fingerprint'].iloc[random_index],filtered_df['fingerprint'].iloc[index])
		filtered_df['jaccard_similarity'].iloc[index] = similarity
		filtered_df['structure_bin_index'].iloc[index] = bin_index

	filtered_df.sort_values(by='jaccard_similarity',inplace=True,ascending=False)
	filtered_df = filtered_df.head(block_size)
	print(filtered_df)
	df.update(filtered_df)
	
print(df)

observed_list = []
frame_list = []

for i in range(0,block_size):
	output_bin_index = 0
	output_string = ""
	counter_list = []
	for j in range(0,number_references):
		bin_identifier = df['structure_bin_index'].iloc[output_bin_index]
		counter_list.append(int(bin_identifier))
		output_bin_index = output_bin_index + i
		output_string = output_string + "," + str(bin_identifier)
	#print(i,output_string)
	bin_frequencies = Counter(counter_list)
	variance = 0
	observed = 0
	for key, value in bin_frequencies.items():
		variance = variance + ((value/number_references - 1/number_references)**2)/number_references
		observed = variance / ( (1/number_references*(1-1/number_references))/number_references * ((len(df)-number_references)/(len(df)-1))    )
	frame_list.append(i)
	observed_list.append(observed)	
	print(variance,observed)
# trendline="rolling", trendline_options=dict(window=5)

fig_4 = px.scatter(x=frame_list,y=observed_list,log_y=True,log_x=True)
fig_4.update_traces(mode='lines')
fig_4.show()
fig_4.write_html("benzene_pocket_convergence_analysis.html")

#################################    Testing end


'''



#print(df['fingerprint'].head())
fingerprint_list = []
fingerprint_list = df['fingerprint'].astype(str)

#tfidf_vectorizer = TfidfVectorizer(max_features=5000)
tfidf_vectorizer = TfidfVectorizer()
X_tfidf = tfidf_vectorizer.fit_transform(fingerprint_list)

print(X_tfidf)

jaccard_sim = jaccard_similarity(df['fingerprint'].iloc[0], df['fingerprint'].iloc[1])
print(f"Jaccard Similarity: {jaccard_sim}")

#cosine_sim = cosine_similarity(X_tfidf[0:1], X_tfidf)
#print(f"Cosine Similarity: {cosine_sim[0][1]}")

#text_features = CountVectorizer(analyzer="word").fit_transform(fingerprint_list)
#print(text_features)

clusterer = KMeans(n_clusters=5)
clusters = clusterer.fit_predict(X_tfidf)
X_pca = PCA(n_components=2).fit_transform(X_tfidf.toarray())

#pca_x_list = []
#pca_y_list = []


with open("pca_analysis_of_pockets.csv","w") as pca_out:
	pca_out.write(f"frame,pca_x,pca_y,cluster\n")
	for i in range(0,len(X_pca)):
#		print(f"{X_pca[i][0]},{X_pca[i][1]},{clusters[i]}")
		pca_out.write(f"{i},{X_pca[i][0]},{X_pca[i][1]},{clusters[i]}\n")
		df.loc[i, 'cluster'] = clusters[i]
		df.loc[i, 'pca_x'] = X_pca[i][0]
		df.loc[i, 'pca_y'] = X_pca[i][1]
		#pca_x_list.append(X_pca[i][0])
		#pca_y_list.append(X_pca[i][1])

df["cluster"] = df["cluster"].astype(str)

fig = px.scatter(df,x='pca_x',y='pca_y',color='frame')
fig.show()
fig.write_html("benzene_pocket_analysis_pca_by_frame.html")

fig_3 = px.scatter(df,x='pca_x',y='pca_y',color='cluster')
fig_3.show()
fig_3.write_html("benzene_pocket_analysis_pca_by_cluster.html")

##############################

xedges = []
yedges = []


for increment in range(-100,100):
        xedges.append(increment*0.01)
        yedges.append(increment*0.01)

H,xbins,ybins = np.histogram2d(df['pca_x'],df['pca_y'], bins=(xedges,yedges))

#print(H)

H_max = np.max(H)

#print(H[0],H_max)


xcenters = (xbins[:-1] + xbins[1:]) / 2
ycenters = (ybins[:-1] + ybins[1:]) / 2

df_2 = pd.DataFrame()
df_2_index = 0

x_list = []
y_list = []

for i in range(0,len(H)):
	for j in range(0,len(H)):
		#print(H[i][j])
		if H[i][j] != 0:
			pmf = -8.314 * 1/4.184 * 1/1000  * 300 * np.log(int(H[i][j])/int(H_max))
			df_2.loc[df_2_index,'xcenter'] = xcenters[i]
			df_2.loc[df_2_index,'ycenter'] = ycenters[j]
			df_2.loc[df_2_index,'pmf'] = pmf
			x_list.append(xcenters[i])
			y_list.append(ycenters[j])
			#print(f"{xcenters[i]},{ycenters[j]},{pmf}")
			df_2_index += 1
		else:
			df_2.loc[df_2_index,'xcenter'] = xcenters[i]
			df_2.loc[df_2_index,'ycenter'] = ycenters[j]
			df_2.loc[df_2_index,'pmf'] = -1
			df_2_index += 1
			
max_value = df_2['pmf'].max()

x_max = max(x_list) + 0.05
x_min = min(x_list) - 0.05

y_max = max(y_list) + 0.05
y_min = min(y_list) - 0.05

print(f"x_max {x_max},x_min {x_min}, y_max {y_max},y_min {y_min}")

df_2['pmf'] = df_2['pmf'].replace(-1, max_value)

df_3 = df_2[ (df_2['xcenter'] < x_max) & (df_2['xcenter'] > x_min) & (df_2['ycenter'] < y_max) & (df_2['ycenter'] > y_min)    ]

# Calculate the min x, min y, max x, and max y values, add a buffer to each one, then create a new dataframe filtered on these values  -- Ben Samudio @ 2025MAR19PTZ2024

fig_2 = px.density_contour(df_3, x="xcenter", y="ycenter", z="pmf", histfunc="min", nbinsx=100, nbinsy=100)
fig_2.update_traces(contours_coloring="heatmap", contours_showlabels = False, colorscale="rainbow")
fig_2.show()
fig_2.write_html("benzene_pocket_analysis_free_energy_landscape.html")


'''
