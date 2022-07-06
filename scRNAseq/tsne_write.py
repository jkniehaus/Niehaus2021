#write tsne coordinates using sklearn's TSNE algorithm.
#perplexity/learning rate were determined via iterative testing
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

sc_pca = pd.read_table('SNIShamcsc_top60pcs.txt',index_col=0)



#when you're ready to write
tsne_structure = TSNE(n_components=2, random_state=123, perplexity=65, learning_rate=5100)
sc_tsne = pd.DataFrame(tsne_structure.fit_transform(sc_pca), index=sc_pca.index)
sc_tsne.to_csv('65perp5100lr.txt',sep='\t')

