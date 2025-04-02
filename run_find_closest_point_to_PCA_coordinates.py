import pandas as pd
import math

x_ref = -0.19 # Global minimum
y_ref = -0.05 # Global minimum

#x_ref = 0.19 # "Saddle" point
#y_ref = 0.01 # "Saddle" point

#x_ref = 0.37 # "meta-stable" state
#y_ref = -0.07 # "meta-stable" state

#x_ref = -0.21
#y_ref = -0.09

df = pd.read_csv("/Users/benjaminsamudio/pca_analysis_of_pockets.csv") 

print(df.info())

old_distance = 100000

for index, row in df.iterrows():
	new_distance = math.sqrt((x_ref - row['pca_x'])**2+(y_ref - row['pca_y'])**2)
	if new_distance < old_distance:
		old_distance = new_distance
		frame = row['frame']
		cluster = row['cluster']

print(old_distance,frame,cluster)
