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


def jaccard_similarity(doc1, doc2):
	words_doc1 = set(doc1.split())
	words_doc2 = set(doc2.split())
	intersection = words_doc1.intersection(words_doc2)
	union = words_doc1.union(words_doc2)
	return len(intersection) / len(union)





df = pd.read_csv("/Users/benjaminsamudio/grid_based_fingerprints.csv")

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

clusterer = KMeans(n_clusters=10)
clusters = clusterer.fit_predict(X_tfidf)
X_pca = PCA(n_components=2).fit_transform(X_tfidf.toarray())

#pca_x_list = []
#pca_y_list = []

for i in range(0,len(X_pca)):
#	print(f"{X_pca[i][0]},{X_pca[i][1]},{clusters[i]}")
	df.loc[i, 'cluster'] = clusters[i]
	df.loc[i, 'pca_x'] = X_pca[i][0]
	df.loc[i, 'pca_y'] = X_pca[i][1]
	#pca_x_list.append(X_pca[i][0])
	#pca_y_list.append(X_pca[i][1])

df["cluster"] = df["cluster"].astype(str)

fig = px.scatter(df,x='pca_x',y='pca_y',color='frame')
fig.show()

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

for i in range(0,len(H)):
	for j in range(0,len(H)):
		#print(H[i][j])
		if H[i][j] != 0:
			pmf = -8.314 * 1/4.184 * 1/1000  * 300 * np.log(int(H[i][j])/int(H_max))
			df_2.loc[df_2_index,'xcenter'] = xcenters[i]
			df_2.loc[df_2_index,'ycenter'] = ycenters[j]
			df_2.loc[df_2_index,'pmf'] = pmf
			#print(f"{xcenters[i]},{ycenters[j]},{pmf}")
			df_2_index += 1
		else:
			df_2.loc[df_2_index,'xcenter'] = xcenters[i]
			df_2.loc[df_2_index,'ycenter'] = ycenters[j]
			df_2.loc[df_2_index,'pmf'] = -1
			df_2_index += 1
			
max_value = df_2['pmf'].max()

df_2['pmf'] = df_2['pmf'].replace(-1, max_value)

# Calculate the min x, min y, max x, and max y values, add a buffer to each one, then create a new dataframe filtered on these values  -- Ben 2025MAR19PTZ2024

fig_2 = px.density_contour(df_2, x="xcenter", y="ycenter", z="pmf", histfunc="min", nbinsx=100, nbinsy=100)
fig_2.update_traces(contours_coloring="heatmap", contours_showlabels = False, colorscale="rainbow")
fig_2.show()
