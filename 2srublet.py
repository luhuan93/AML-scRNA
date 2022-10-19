import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

input_dir = '/data/AMLscRNA/1026/outs/filtered_feature_bc_matrix'  
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix) 
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.call_doublets(threshold=0.25)
scrub.plot_histogram()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
print (scrub.detected_doublet_rate_) 0.0
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P1b_doublet.txt', index=False,header=True)
out_df.head()

input_dir = '/data/AMLscRNA/1163/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
print (scrub.detected_doublet_rate_) 0.02999210734017364
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P1a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1059/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
print (scrub.detected_doublet_rate_) 0.021804966686856452
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P2b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1176/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
print (scrub.detected_doublet_rate_) 0.045276690888764674
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P2a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1077/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.00023772732675621063
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P3b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1191/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.032353688489878886
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P3a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1266/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.00018955549237039144
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P4b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1421/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.02527948433880552
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P4a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1152/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.00072992700729927
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P5b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1320/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.022199627181833586
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P5a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1363/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.0005136546528550637
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P6b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1503/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.04450065579913809
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P6a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1309/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.0
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P7b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/1449/outs/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.05650657790458076
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P7a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/YMF-B-1677/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])


scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.018304530115553736
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P8b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/YMF-A-1817/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.005709065996802923
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P8a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/HZW-B-2238/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.0066518847006651885
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P9b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/HZW-A-2381/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.023224653863331845
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P9a_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/ZYM-B-2215/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.024756285274499742
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P10b_doublet.txt', index=False,header=True)

input_dir = '/data/AMLscRNA/ZYM-A-2328/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(input_dir + '/barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.018686482351655556
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/scrublet/P10a_doublet.txt', index=False,header=True)

import pathlib
import glob
input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/26_1')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 9.986019572598362e-05
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/26_1_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/26_2')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.00045202847779410104
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/26_2_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/26_3')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.019801980198019802
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/26_3_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/27_1')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.001155067860236789
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/27_1_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/27_2')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.00617828773168579
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/27_2_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/27_3')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.0008821453775582216
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/27_3_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/28_1')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.002998092123194331
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/28_1_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/28_2')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.015486725663716814
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/28_2_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/28_3')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.008824093612123537
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/28_3_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/29_1')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.0014892032762472078
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/29_1_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/29_2')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.012307692307692308
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/29_2_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/29_3')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.010676873489121675
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/29_3_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/30_1')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.005003001801080648
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/30_1_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/30_2')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.000955718381650207
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/30_2_doublet.txt', index=False,header=True)

input_dir = pathlib.Path('/data/AMLscRNA/BM/blood/30_3')
counts_matrix = scipy.io.mmread(list(input_dir.glob("*.mtx.gz"))[0]).T.tocsc()
genes = np.array(scr.load_genes(list(input_dir.glob("*.tsv"))[0], delimiter='\t', column=1)) #要解压
out_df = pd.read_csv(list(input_dir.glob("*.tsv.gz"))[0], header = None, index_col=None, names=['barcode'])

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.show()
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()
print (scrub.detected_doublet_rate_) 0.01025568980050576
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv('/data/AMLscRNA/BM/blood/scrublet/30_3_doublet.txt', index=False,header=True)
