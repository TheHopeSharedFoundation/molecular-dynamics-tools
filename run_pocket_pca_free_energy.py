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

df = pd.read_csv("/Users/benjaminsamudio/grid_based_fingerprints-TEST.csv")
fingerprint_list = []
fingerprint_list = df['fingerprint'].astype(str)

tfidf_vectorizer = TfidfVectorizer()
X_tfidf = tfidf_vectorizer.fit_transform(fingerprint_list)
array_tfidf = X_tfidf.toarray()

#cosine_sim = cosine_similarity(array_tfidf, array_tfidf)

clusterer = KMeans(n_clusters=5)
clusters = clusterer.fit_predict(X_tfidf)

X_pca = PCA(n_components=2).fit_transform(X_tfidf.toarray())


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
