import loompy
loompy.combine(files=['/data/AMLscRNA/1026/velocyto/1026.loom',
'/data/AMLscRNA/1163/velocyto/1163.loom',
'/data/AMLscRNA/1059/velocyto/1059.loom',
'/data/AMLscRNA/1176/velocyto/1176.loom',
'/data/AMLscRNA/1077/velocyto/1077.loom',
'/data/AMLscRNA/1191/velocyto/1191.loom',
'/data/AMLscRNA/1266/velocyto/1266.loom',
'/data/AMLscRNA/1421/velocyto/1421.loom',
'/data/AMLscRNA/1152/velocyto/1152.loom',
'/data/AMLscRNA/1320/velocyto/1320.loom',
'/data/AMLscRNA/1363/velocyto/1363.loom',
'/data/AMLscRNA/1503/velocyto/1503.loom',
'/data/AMLscRNA/1309/velocyto/1309.loom',
'/data/AMLscRNA/1449/velocyto/1449.loom',
'/data/AMLscRNA/HZW-A-2381/velocyto/HZW-A-2381.loom',
'/data/AMLscRNA/HZW-B-2238/velocyto/HZW-B-2238.loom',
'/data/AMLscRNA/YMF-A-1817/velocyto/YMF-A-1817.loom',
'/data/AMLscRNA/YMF-B-1677/velocyto/YMF-B-1677.loom',
'/data/AMLscRNA/ZYM-A-2328/velocyto/ZYM-A-2328.loom',
'/data/AMLscRNA/ZYM-B-2215/velocyto/ZYM-B-2215.loom'
], key="Accession",output_file="/data/AMLscRNA/mergeblood/merge.loom")
tot=sc.read_loom('/data/AMLscRNA/mergeblood/merge.loom', sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
tot.obs.index = [i.split('-')[0] for i in tot.obs.index]


import os
import scvelo as scv
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 
os.chdir('/data/AMLscRNA/mergeblood/20220501/undiff')
tot=scv.read_loom('/data/AMLscRNA/mergeblood/merge.loom', sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
adata = scv.read("/data/AMLscRNA/mergeblood/20220501/undiff/combined_ref_10p_0501.h5ad", cache=True)
undiff =scv.read("/data/AMLscRNA/mergeblood/20220501/undiff/undiff0501.h5ad", cache=True)
adata= adata[adata.obs['Cohort']=="1",:]
adata.obs.index = [i.split('-')[0] for i in adata.obs.index]
undiff.obs.index = [i.split('-')[0] for i in undiff.obs.index]
undiff.obs['celltype']=adata.obs['predictions']
undiff.obsm['X_umap']=adata.obsm['X_umap']

import re
r = re.compile('P[1-6]$')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
response=undiff[vmatch(undiff.obs['sample']),:]
response.obs['sample'].unique()

r = re.compile('P[7-9]$|P10')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
unres=undiff[vmatch(undiff.obs['sample']),:]
unres.obs['sample'].unique()
RESa=response[response.obs['state']== 'after']
RESb=response[response.obs['state']== 'before']
UNa=unres[unres.obs['state']== 'after']
UNb=unres[unres.obs['state']== 'before']

RESadata = scv.utils.merge(RESa, tot)
scv.pp.filter_and_normalize(RESadata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(RESadata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(RESadata, n_jobs=32)
scv.tl.velocity(RESadata, mode='dynamical')
scv.tl.velocity_graph(RESadata)
scv.tl.paga(RESadata, groups='celltype')
sns.set(rc={"figure.figsize":(6, 4)})
scv.pl.paga(RESadata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=2, arrowsize = 15,save='RESa_scvelo_PAGA1.pdf')

RESbdata = scv.utils.merge(RESb, tot)
scv.pp.filter_and_normalize(RESbdata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(RESbdata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(RESbdata, n_jobs=32)
scv.tl.velocity(RESbdata, mode='dynamical')
scv.tl.velocity_graph(RESbdata)
scv.tl.paga(RESbdata, groups='celltype')
sns.set(rc={"figure.figsize":(6, 4)})
scv.pl.paga(RESbdata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=2, arrowsize = 15,save='RESb_scvelo_PAGA1.pdf')

UNbdata = scv.utils.merge(UNb, tot)
scv.pp.filter_and_normalize(UNbdata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(UNbdata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(UNbdata, n_jobs=32)
scv.tl.velocity(UNbdata, mode='dynamical')
scv.tl.velocity_graph(UNbdata)
scv.tl.paga(UNbdata, groups='celltype')
sns.set(rc={"figure.figsize":(6, 4)})
scv.pl.paga(UNbdata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=2, arrowsize = 15,save='UNb_scvelo_PAGA1.pdf')

UNadata = scv.utils.merge(UNa, tot)
scv.pp.filter_and_normalize(UNadata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(UNadata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(UNadata, n_jobs=32)
scv.tl.velocity(UNadata, mode='dynamical')
scv.tl.velocity_graph(UNadata)
scv.tl.paga(UNadata, groups='celltype')
sns.set(rc={"figure.figsize":(6, 4)})
scv.pl.paga(UNadata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=2, arrowsize = 15,save='UNa_scvelo_PAGA1.pdf')

