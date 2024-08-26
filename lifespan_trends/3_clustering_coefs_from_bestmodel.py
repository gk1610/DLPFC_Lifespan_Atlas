
import pandas as pd
import pyarrow.feather as feather
import pandas as pd
import patsy
from sklearn.cluster import KMeans
from sklearn import metrics
import numpy as np

df = feather.read_feather('lifespan_trends/best_model/best_model_coefs.feather')
df.set_index('celltype_ID', inplace=True)
X_no_intercept = df

kmeans = KMeans(n_clusters=10, random_state=42,max_iter=1000)
clusters = kmeans.fit_predict(X_no_intercept)
df['cluster'] = clusters
df.to_csv('lifespan_trends/best_model/clustered_data_nclust10.csv', index=True)
