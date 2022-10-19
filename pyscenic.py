import os
os.chdir('/data/AMLscRNA/mergeblood/20220501/pyscenic/raw')
#pyscenic
import loompy as lp
import pandas as pd
AMLy = sc.read_loom("/data/AMLscRNA/AML_10pC.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
import scanpy as sc
AML=sc.read_h5ad('/data/AMLscRNA/mergeblood/20220501/AML10_with_raw.h5ad')
row_attrs = {
"Gene": np.array(AMLy.var_names),
}
col_attrs = {
"CellID": np.array(AMLy.obs_names),
"nGene": np.array(np.sum(AMLy.X.transpose()>0, axis=0)).flatten(),
"nUMI": np.array(np.sum(AMLy.X.transpose(),axis=0)).flatten(),
}
lp.create("AML_pyscenic_raw.loom",AMLy.X.transpose(),row_attrs,col_attrs)

import scanpy as sc
AML=sc.read_h5ad('AML10_with_raw.h5ad')

'''
/share/home/luhuan//miniconda3/bin/arboreto_with_multiprocessing.py \
-m grnboost2 \
-o ./adj.AML_raw.tsv \
--num_workers 90 \
-t ./AML_pyscenic_raw.loom \
./hs_hgnc_tfs.txt
s
/share/home/luhuan//miniconda3/bin/pyscenic ctx \
adj.AML_raw.tsv \
../hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname ../motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ./AML_pyscenic_raw.loom \
--mode "dask_multiprocessing" \
--output ./AML_reg_raw.csv \
--num_workers 100 \
--mask_dropouts

/share/home/luhuan//miniconda3/bin/pyscenic aucell \
./AML_pyscenic_raw.loom \
./AML_reg_raw.csv \
--output ./AML_SCENIC_raw.loom \
--num_workers 100
'''

#create neutro_SCENIC_integrate.loom
lf = lp.connect("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/AML_SCENIC_raw.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
#lf.close()
import umap
# UMAP
runUmap = umap.UMAP(n_neighbors=20, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_umap.txt", sep='\t')

meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'scenic_umap.txt', sep='\t', header=0, index_col=0 )

auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
        pd.DataFrame(AML.obsm['X_umap'],index=AML.obs.index)[0] ,
        pd.DataFrame(AML.obsm['X_pca'],index=AML.obs.index)[0] ,
        dr_umap['X']
    ], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
        pd.DataFrame(AML.obsm['X_umap'],index=AML.obs.index)[1] ,
        pd.DataFrame(AML.obsm['X_pca'],index=AML.obs.index)[1] ,
        dr_umap['Y']
    ], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3']

metaJson = {}

metaJson['embeddings'] = [
    {
        "id": 1,
        "name": f"Scanpy UMAP  (highly variable genes)"
    },
    {
        "id": 2,
        "name": "Scanpy PC1/PC2"
    },
    {
        "id": 3,
        "name": "SCENIC AUC UMAP"
    },
]

metaJson["clusterings"] = [{
            "id": 0,
            "group": "Scanpy",
            "name": "Scanpy leiden 0.89",
            "clusters": [],
        }]

metaJson["metrics"] = [
        {
            "name": "nUMI"
        }, {
            "name": "nGene"
        }, {
            "name": "Percent_mito"
        }
]
import numpy as np
metaJson["annotations"] = [
    {
        "name": "Ledi_clusters_Scanpy",
        "values": list(set( AML.obs['leiden_089'].astype(np.str) ))
    },]    	
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in AML.obs['leiden_089']])) + 1):
    clustDict = {}
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i + 1}'
    metaJson['clusterings'][0]['clusters'].append(clustDict)
    
clusterings = pd.DataFrame()
clusterings["0"] = AML.obs['leiden_089'].values.astype(np.int64)

def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr
    
col_attrs = {
    "CellID": np.array(AML.obs.index),
    "nUMI": np.array(AML.obs['n_counts'].values),
    "nGene": np.array(AML.obs['nFeature_RNA'].values),
    "leiden_clusters_Scanpy": np.array( AML.obs['leiden_089'].values ),
    "Percent_mito": np.array(AML.obs['percent.mt'].values),
    "Embeddings_X": dfToNamedMatrix(Embeddings_X),
    "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
    "RegulonsAUC": dfToNamedMatrix(auc_mtx),
    "Clusterings": dfToNamedMatrix(clusterings),
    "celltype": np.array(AML.obs['CellTypeN'].values)
}

row_attrs = {
    "Gene": lf.ra.Gene,
    "Regulons": regulons,
}

attrs = {
    "title": "sampleTitle",
    "MetaData": json.dumps(metaJson),
    "Genome": 'hg38',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')
lp.create(
    filename = "/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/AML_SCENIC_integrate.loom" ,
    layers=lf[:,:],
    row_attrs=row_attrs, 
    col_attrs=col_attrs, 
    file_attrs=attrs
)
lf.close()   

import json
import zlib
import base64
import loompy as lp
import pandas as pd
lf = lp.connect("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/AML_SCENIC_integrate.loom", mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
dr_names = [
    meta['embeddings'][0]['name'].replace(" ","_")
]

drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )

dr=[]
for i in range( len(drx.columns) ):
    dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
   
# rename columns:
for i,x in enumerate( dr ):
    x.columns = ['X','Y']
    
def colorMap( x, palette=my_palette ):
    import natsort
    from collections import OrderedDict
    #
    n=len(set(x))
    cpalette = sns.color_palette(palette,n_colors=n )
    cdict = dict( zip( list(set(x)), cpalette ))
    cmap = [ cdict[i] for i in x ]
    cdict = OrderedDict( natsort.natsorted(cdict.items()) )
    return cmap,cdict

def drplot( dr, colorlab, ax, palette=my_palette, title=None):
    cmap,cdict = colorMap( colorlab, palette )
    for lab,col in cdict.items():  
        ix = colorlab.loc[colorlab==lab].index
        ax.scatter( dr['X'][ix], dr['Y'][ix], c=[col]*len(ix), alpha=0.7, s=0.8,label=lab, edgecolors='none')
    if( title is not None ):
        ax.set_title(title, fontsize='xx-large');
    #
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)    

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10})
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(23,10), dpi=150 )
cellAnnot = pd.concat(
    [
        pd.DataFrame( lf.ca.leiden_clusters_Scanpy, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.celltype, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.Percent_mito, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
    ],
    axis=1
)
cellAnnot.columns = [
 'ClusterID',
 'celltype',
 'Percent_mito',
 'nGene',
 'nUMI']
import seaborn as sns
my_palette={'B':'#023fa5', 'CD14+mono':'#7d87b9', 'CD16+mono':'#bec1d4', 'Early_ery':'#d6bcc0', 'Ery':'#bb7784', 'GMP':'#8e063b', 'HSPC':'#4a6fe3', 'M2_Macro':'#8595e1', 'MEP':'#b5bbe3', 'NK':'#e6afb9', 'NKT':'#e07b91', 'Plasma':'#d33f6a', 'Platelets':'#11c638', 'Pre-B':'#8dd593', 'T':'#c6dec7', 'cDC':'#ead3c6', 'pDC':'#f0b98d', 'Stromal':'#ef9708','LSPC_Cycle':'#FF7F0E', 'LSPC_Primed':'#279E68', 'LSPC_Quiescent':'#D62728','Mono_like':'#AA40FC', 'ProMono_like':'#8C564B', 'cDC_like':'#E377C2', 'GMP_like':'#1F77B4'}
drplot( dr[0], colorlab=cellAnnot['celltype'], ax=ax1, palette=colormap, title='Highly variable genes - UMAP' )
drplot( dr[2], colorlab=cellAnnot['celltype'], ax=ax2, palette=colormap, title='SCENIC AUC - UMAP' )
ax2.legend( bbox_to_anchor=(1,1), ncol=1,  markerscale=10,title_fontsize='xx-large',fontsize='xx-large', frameon=False, title="celltype")
plt.tight_layout()
plt.savefig("umap-hvg-scenic-louvain.png", dpi=600, bbox_inches = "tight")   


import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import scanpy as sc
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from adjustText import adjust_text
from pyscenic.binarization import binarize
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
sig = load_signatures('/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/AML_reg_raw.csv')
adata = add_scenic_metadata(AML, auc_mtx, sig)

adata.obs['celltype1']=[i +'_'+ j for i, j in zip(adata.obs['type'],adata.obs['CellTypeN'])]
df_obs = adata.obs
order=['LSPC_Quiescent','LSPC_Primed',  'LSPC_Cycle', 'GMP_like', 'ProMono_like', 'Mono_like', 'cDC_like']
order1=['RESb','RESa','UNb', 'UNa']
order2=[a +'_'+b for a in order1 for b in order ]
topreg = []
for i,c in enumerate(order2):
    topreg.extend(
        list(selec.T[c].sort_values(ascending=False)[:10].index)
    )
topreg=['Regulon('+i + ')' for i in topreg]
df_obs = df_obs[df_obs['celltype1'].isin(order2) ]
signature_column_names= list(set(topreg))
df_scores = df_obs[signature_column_names + ['celltype1']]
df_results = ((df_scores.groupby(by='celltype1').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results=pd.pivot_table(data=df_results,index='celltype1', columns='regulon', values='Z')
fig, ax1 = plt.subplots(1, 1, figsize=(20, 10))
ax=sns.clustermap(df_results,  linewidths=.7, cbar=True, figsize=(20, 10),square=True, linecolor='black', cmap="RdBu_r")
plt.savefig("split_tumor_regulontop10.pdf")


df_obs = adata.obs[adata.obs['type']=="RESb"]
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['CellTypeN']]
df_results = ((df_scores.groupby(by='CellTypeN').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results=pd.pivot_table(data=df_results,index='CellTypeN', columns='regulon', values='Z')
k=rss_celltypeRESb
k=rss_celltypeRESa
k=rss_celltypeUNb
k=rss_celltypeUNa
order=['LSPC_Quiescent','LSPC_Primed',  'LSPC_Cycle', 'GMP_like', 'ProMono_like', 'Mono_like', 'cDC_like']
selec=k.loc[order]
topreg = []
for i,c in enumerate(order):
    topreg.extend(
        list(selec.T[c].sort_values(ascending=False)[:10].index)
    )
RESb_regulontop10 =df_results[set(topreg)]

regulon_filter=pd.read_csv("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/regulon_filter_tumor", delimiter='\t')
regulon_filter=pd.pivot_table(data=regulon_filter,index='celltype', columns='regulon', values='Z')
regulon_filter=regulon_filter.reindex(index= order)
fig, ax1 = plt.subplots(1, 1, figsize=(20, 10))
sns.heatmap(regulon_filter, ax=ax1, linewidths=.7, cbar=True, square=True, linecolor='black', 
            cmap="RdBu_r")
ax=sns.clustermap(regulon_filter,  linewidths=.7, cbar=True, figsize=(20, 5),square=True, linecolor='black', cmap="RdBu_r")
plt.savefig("cellType-regulon-zscore_tumor.pdf")
plt.show()

cell.prop2<-read.table("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/regulon0525.csv",sep="\t",header = T)

cell=pd.read_csv("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/regulon0525.csv", delimiter='\t')
regulon_filter=pd.pivot_table(data=cell,index='celltype', columns='regulon', values='Z')
regulon_filter.to_csv('regulon_matrix0525.csv', index=True,header=True,sep="\t")
'''
cell.prop2$celltype = factor(cell.prop2$celltype, levels=c('LSPC_Quiescent','LSPC_Primed',  'LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like'))
ggplot(cell.prop2, aes( regulon,celltype)) +
    geom_tile(aes(fill = z_score, colour = "white")) +
    scale_fill_gradient2(low = "blue", high = "red",mid="white",midpoint = 0)+scale_x_continuous(limits = c(0, 10)
'''

df_results[(df_results.Z >= 2.0)].sort_values('Z', ascending=False).head()
df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 2].sort_values('Z', ascending=False),
                           index='CellTypeN', columns='regulon', values='Z')
fig, ax1 = plt.subplots(1, 1, figsize=(20, 16))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray', 
            cmap="YlGnBu", annot_kws={"size": 2})
ax1.set_ylabel('')
plt.savefig("cellType-regulon-zscore.pdf", dpi=600, bbox_inches = "tight")
plt.show()


typem=AML10.obs['type']

temauc=pd.merge(auc_mtx,typem,how='left', left_index=True, right_index=True)
cellAnnot['celltype1']=[i +'_'+ j for i, j in zip(AML10.obs['type'],AML10.obs['CellTypeN'])]

auc_mtxRESb=temauc[temauc['type']=="RESb"].drop(['type'],axis=1)
cellAnnotRESb=AML10.obs[AML10.obs['type']=="RESb"]['CellTypeN']
rss_celltypeRESb = regulon_specificity_scores( auc_mtxRESb, cellAnnotRESb)

auc_mtxRESa=temauc[temauc['type']=="RESa"].drop(['type'],axis=1)
cellAnnotRESa=AML10.obs[AML10.obs['type']=="RESa"]['CellTypeN']
rss_celltypeRESa = regulon_specificity_scores( auc_mtxRESa, cellAnnotRESa)

auc_mtxUNb=temauc[temauc['type']=="UNb"].drop(['type'],axis=1)
cellAnnotUNb=AML10.obs[AML10.obs['type']=="UNb"]['CellTypeN']
rss_celltypeUNb = regulon_specificity_scores( auc_mtxUNb, cellAnnotUNb)

auc_mtxUNa=temauc[temauc['type']=="UNa"].drop(['type'],axis=1)
cellAnnotUNa=AML10.obs[AML10.obs['type']=="UNa"]['CellTypeN']
rss_celltypeUNa = regulon_specificity_scores( auc_mtxUNa, cellAnnotUNa)

rss_celltype = regulon_specificity_scores( auc_mtx, cellAnnot['celltype'] )
rss_celltypeall = regulon_specificity_scores( auc_mtx, cellAnnot['celltype1'] )
cats = sorted( list(set(cellAnnot['celltype'])) )
cats = sorted( list(set(cellAnnotRESb)) )
cats = sorted( list(set(cellAnnotRESa)) )
cats = sorted( list(set(cellAnnotUNb)) )
cats = sorted( list(set(cellAnnotUNa)) )
groupby_order=['HSPC',  'GMP','MEP','Pre-B','B','Plasma', 'T', 'NK',  'CD14+mono', 'CD16+mono','Platelets', 'pDC',  'cDC','Early_ery', 'Ery', 'Stroma', 'LSPC_Primed', 'LSPC_Quiescent', 'LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
C = {val:ix for ix,val in enumerate(groupby_order)}	
sort1=sorted(cats, key=groupby_order.index)	

k=rss_celltypeRESb
k=rss_celltypeRESa
k=rss_celltypeUNb
k=rss_celltypeUNa
fig = plt.figure(figsize=(15, 12))
for c,num in zip(sort1, range(1,len(sort1)+1)):
    x=k.T[c]
    ax = fig.add_subplot(5,5,num)
    plot_rss(k, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large',
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig("RSS-top5.pdf", dpi=150, bbox_inches = "tight")
plt.show()


import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display
from collections import OrderedDict

#选择5个进行展示
topreg = []
for i,c in enumerate(sort):
    topreg.extend(
        list(rss_celltype.T[c].sort_values(ascending=False)[:10].index)
    )

topreg = list(set(topreg))
topreg = pd.DataFrame(topreg)
topreg.to_csv( "top5_reg.txt", sep='\t',index=False)
cell=pd.read_csv("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/top5_reg.bed", delimiter='\t')
regulon_filter=pd.pivot_table(data=cell,index='celltype', columns='regulon', values='Z')
regulon_filter.to_csv('regulon_matrix_top5.csv', index=True,header=True,sep="\t")

regu = AML.obs.iloc[:, 25:]
regu.to_csv('regulon_matrix_RAS_raw.csv', index=False,header=False,sep="\t")
pd.DataFrame(age_sex.columns).to_csv('regulon_matrix_RAS_raw_column.csv', index=False,header=False,sep="\t")
pd.DataFrame(AML.obs[['sam','state','CellTypeN']]).to_csv('regulon_matrix_RAS_raw_celltype.csv', index=False,header=False,sep="\t")
'''
#diff activity regulon
regulonM<-read.table("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/regulon_matrix_RAS_raw.csv",head=T,sep="\t")
regulonM<-as.matrix(regulonM)
fData<-read.table("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/regulon_matrix_RAS_raw_column.csv",head=T,sep="\t")
fData<-data.frame(fData)
cData<-read.table("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/regulon_matrix_RAS_raw_celltype.bed",head=T,sep="\t")
cData<-data.frame(cData)
sca <- FromMatrix(t(regulonM), cData, fData)
filterCrit <- with(colData(sca), celltype=="LSPC_Quiescent"& type =="UN")
tem <- subset(sca,filterCrit)
cond<-factor(colData(tem)$state)
cond<-relevel(cond,"before")
colData(tem)$state<-cond
zlmCond <- zlm(~state , tem)
summaryCond <- summary(zlmCond, doLRT='stateafter') 
print(summaryCond, n=4)
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='stateafter' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='stateafter' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05], as.data.table(mcols(tem)), by='primerid')
setorder(fcHurdleSig, fdr)
write.table(fcHurdleSig,"LSPC_Quiescent_regulon_diff",sep = "\t",quote=F,row.names = F)
'''
AML.obs['type']=AML.obs['type'].cat.add_categories(["RES","UN"])
import re
r = re.compile('P[1-6]')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
AML.obs.loc[vmatch(AML.obs['orig.ident']),'type']="RES"
r = re.compile('P[7-9]|P10')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
AML.obs.loc[vmatch(AML.obs['orig.ident']),'type']="UN"

sns.set(rc={"figure.figsize":(8, 5)}) 
sns.set(style="ticks")

order=['RESb','RESa','UNb','UNa']
tem = AML[(AML.obs['CellTypeN'] == 'HSPC'),:].obs[['type','Regulon(NFKB2_(+))']]
tem.to_csv('HSPC_NFKB2_RAS.csv', index=False,header=False,sep="\t")
tem = AML[(AML.obs['CellTypeN'] == 'HSPC'),:].obs[['type','Regulon(ATF3_(+))']]
tem = AML[(AML.obs['CellTypeN'] == 'HSPC'),:].obs[['type','Regulon(BCL3_(+))']]
tem.to_csv('HSPC_NFKB2_RAS.csv', index=False,header=False,sep="\t")
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed'),:].obs[['type','Regulon(TP53_(+))']]
tem.to_csv('LSPC_Primed_TP53_RAS.csv', index=False,header=False,sep="\t")
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed')|(AML.obs['CellTypeN'] == 'LSPC_Cycle')|(AML.obs['CellTypeN'] == 'LSPC_Quiescent')|(AML.obs['CellTypeN'] == 'GMP_like'),:].obs[['type','Regulon(HOXA7_(+))']]
tem.to_csv('4tumor_HOXA7_RAS.csv', index=False,header=False,sep="\t")
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed')|(AML.obs['CellTypeN'] == 'LSPC_Cycle')|(AML.obs['CellTypeN'] == 'LSPC_Quiescent')|(AML.obs['CellTypeN'] == 'GMP_like'),:].obs[['type','Regulon(RFX8_(+))']]
tem.to_csv('4tumor_RFX8_RAS.csv', index=False,header=False,sep="\t")
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed')|(AML.obs['CellTypeN'] == 'LSPC_Cycle')|(AML.obs['CellTypeN'] == 'LSPC_Quiescent')|(AML.obs['CellTypeN'] == 'GMP_like'),:].obs[['type','Regulon(BHLHE40_(+))']]
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed')|(AML.obs['CellTypeN'] == 'LSPC_Cycle')|(AML.obs['CellTypeN'] == 'LSPC_Quiescent')|(AML.obs['CellTypeN'] == 'GMP_like'),:].obs[['type','Regulon(SMAD1_(+))']]
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed')|(AML.obs['CellTypeN'] == 'LSPC_Cycle')|(AML.obs['CellTypeN'] == 'LSPC_Quiescent')|(AML.obs['CellTypeN'] == 'GMP_like'),:].obs[['type','Regulon(JUNB_(+))']]
tem = AML[(AML.obs['CellTypeN'] == 'LSPC_Primed')|(AML.obs['CellTypeN'] == 'LSPC_Cycle')|(AML.obs['CellTypeN'] == 'LSPC_Quiescent')|(AML.obs['CellTypeN'] == 'GMP_like'),:].obs[['type','Regulon(KLF8_(+))']]

fig=sns.violinplot(x='type', y='Regulon(NFKB2_(+))', data=tem,order=order,inner="quartile",save="NFKB2_boxplot.pdf")
fig.set(title="HSPC")
plt.savefig("NFKB2_boxplot.pdf")
fig=sns.violinplot(x='type', y='Regulon(BCL3_(+))', data=tem,order=order,inner="quartile",save="BCL3_boxplot.pdf")
fig.set(title="HSPC")
plt.savefig("BCL3_boxplot.pdf")
fig=sns.violinplot(x='type', y='Regulon(ATF3_(+))', data=tem,order=order,inner="quartile",save="ATF3_boxplot.pdf")
fig.set(title="HSPC")
plt.savefig("ATF3_boxplot.pdf")
plt.show()
fig=sns.violinplot(x='type', y='Regulon(RFX8_(+))', data=tem,order=order,inner="quartile",cut=0)
fig.set(title="4 tumor type")
plt.savefig("RFX8_violinplot.pdf")
plt.show()

fig=sns.violinplot(x='type', y='Regulon(HOXA7_(+))', data=tem,order=order,inner=None,cut=0)
fig.set(title="4 tumor type")
plt.savefig("HOXA7_violinplot.pdf")
plt.show()

fig=sns.violinplot(x='type', y='Regulon(BHLHE40_(+))', data=tem,order=order,inner=None,cut=0)
fig.set(title="4 tumor type")
plt.savefig("BHLHE40_violinplot.pdf")
plt.show()

fig=sns.violinplot(x='type', y='Regulon(SMAD1_(+))', data=tem,order=order,inner=None,cut=0)
fig.set(title="4 tumor type")
plt.savefig("SMAD1_violinplot.pdf")
plt.show()

fig=sns.violinplot(x='type', y='Regulon(JUNB_(+))', data=tem,order=order,inner=None,cut=0)
fig.set(title="4 tumor type")
plt.savefig("JUNB_violinplot.pdf")
plt.show()

fig=sns.violinplot(x='type', y='Regulon(KLF8_(+))', data=tem,order=order,inner=None,cut=0)
fig.set(title="4 tumor type")
plt.savefig("KLF8_violinplot.pdf")
plt.show()

y="Regulon(BHLHE40_(+))"
y="Regulon(SMAD1_(+))"
y="Regulon(JUNB_(+))"
y="Regulon(KLF8_(+))"
tem=AML[(AML.obs['celltype'].str.startswith("undiff"))&(AML.obs['state'] == 'before'),:].obs[['type',y, 'CellTypeN']]
sort1=sorted(tem['CellTypeN'].unique().tolist(), key=C.get) 
import seaborn as sns
sns.set(rc={"figure.figsize":(16, 16)})
sns.set(style="ticks")
order = ['RESb','UNb']
g=sns.catplot(x = "type", y = y, hue = "type",linewidth=0.4, cut=0,kind = 'violin',inner=None,order=order, col="CellTypeN", data = tem,col_order=sort1,height=3,aspect=0.8,col_wrap=7,sharey=True,dodge=False)
g.set_xticklabels(rotation=45, fontsize = 12)
for (row_key),ax in g.axes_dict.items():
    ax.set_title(f"{row_key}", rotation=0, fontsize = 12)
    ax.set_xlabel(xlabel = None)
    ax.set_xticklabels(labels=order,rotation=45)
       
ax.set_xticklabels(order, fontsize=12) 
g.tight_layout(pad=1, w_pad=1, h_pad=1) 
g.fig.subplots_adjust(wspace=0)
g.fig.subplots_adjust(hspace=0)
g.despine(top=False,right=False,left=False,bottom=False)
import matplotlib.pyplot as plt
plt.savefig("celltype_"+y+"_violin.pdf")


adjacencies = pd.read_csv("/data/AMLscRNA/mergeblood/20220501/pyscenic/raw/adj.AML_raw.tsv", index_col=False, sep='\t')
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.type).T
from pyscenic.utils import modules_from_adjacencies
modules = list(modules_from_adjacencies(adjacencies, exprMat))
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)
tf = 'HOXA7'
tf = 'RFX8'
tf_mods = [ x for x in modules if x.transcription_factor==tf ]

for i,mod in enumerate( tf_mods ):
    print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
print( f'{tf} regulon: {len(regulons[tf+"_(+)"])} genes' )

for i,mod in enumerate( tf_mods ):
    with open( tf+'_module_'+str(i)+'.txt', 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)
            
with open( tf+'_regulon.txt', 'w') as f:
    for item in regulons[tf+'_(+)']:
        f.write("%s\n" % item)



















































import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scanpy.plotting._tools.scatterplots import plot_scatter 新版本没有
import seaborn as sns



meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)
    

 

    
lf.close()
	
for i in range(max(set([int(x) for x in AML.obs['CellTypeN']])) + 1)