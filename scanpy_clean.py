import numpy as np
import pandas as pd
import scanpy as sc
import sc_toolbox.api as sct
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import scanpy.external as sce
import matplotlib.pyplot as plt
import os
os.chdir('/data/AMLscRNA/mergeblood/20220501')
AMLy = sc.read_loom("/data/AMLscRNA/AML_10pC.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
heal = sc.read_loom("/data/AMLscRNA/BM/blood/blood5n.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
AML = AMLy.concatenate(heal,batch_key='batch') 
AML.obs['batch']=AML.obs['batch'].cat.add_categories(['health','patient']) 
AML.obs.loc[AML.obs['batch']=='0','batch']='patient'
AML.obs.loc[AML.obs['batch']=='1','batch']='health'
AML.obs['batch']=AML.obs['batch'].cat.remove_unused_categories() 

pd.crosstab(columns=AML.obs['batch'],index =AML.obs['orig.ident']).to_csv('/data/AMLscRNA/mergeblood/sample_check0405.csv', index=True,header=True,sep="\t")
# ribosomal genes
AML.var['ribo'] = AML.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
AML.var['hb'] = AML.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(AML, qc_vars=['ribo','hb'], percent_top=None, log1p=False, inplace=True)
import seaborn as sns
sns.set(rc={"figure.figsize":(8, 4)}) 
sns.set(style="ticks")
sc.pl.violin(AML, ['n_genes_by_counts'],jitter=0.4, groupby = 'orig.ident', rotation= 45,size=0,figsize=(12,4),save='n_genes_by_counts.png')
sc.pl.violin(AML, ['total_counts'],jitter=0.4, groupby = 'orig.ident', rotation= 45,size=0,figsize=(12,4),save='total_counts.png')
sc.pl.violin(AML, ['percent.mt'],jitter=0.4, groupby = 'orig.ident', rotation= 45,size=0.1,save='percent.mt.png')
sc.pl.violin(AML, ['pct_counts_ribo'],jitter=0.4, groupby = 'orig.ident', rotation= 45,size=0.1,save='pct_counts_ribo.png')
sc.pl.violin(AML, ['pct_counts_hb'],jitter=0.4, groupby = 'orig.ident', rotation= 45,size=0.1,save='pct_counts_hb.png')

sc.pl.highest_expr_genes(AML, n_top=20,save='Top20_dedoublet.png')

mito_genes = AML.var_names.str.startswith('MT-')
AML.obs['percent_mt2'] = np.sum(
    AML[:, mito_genes].X, axis=1).A1 / np.sum(AML.X, axis=1).A1
AML.obs['n_counts'] = AML.X.sum(axis=1).A1

sc.pl.scatter(AML, x='total_counts', y='percent.mt', color="orig.ident",save='count_MT_scatter_dedoublet.png')
  	
cell_cycle_genes = [x.strip() for x in open('/data/AMLscRNA/cell_score')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in AML.var_names]

sc.pp.normalize_total(AML, target_sum=1e4)
sc.pp.log1p(AML)
AML.raw = AML
sc.tl.score_genes_cell_cycle(AML, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(AML, ['S_score', 'G2M_score'],jitter=0.4, groupby = 'orig.ident', rotation=45,size=0.1,save='cellcyle_dedoublet.png')
AML.write_h5ad('AML10pC_blood5n.h5ad')
AML = sc.read_h5ad('/data/AMLscRNA/mergeblood/AML10pC_blood5n.h5ad') 
AML2 = AML.raw.to_adata() 
sc.pp.highly_variable_genes(AML2, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key="batch",subset=True)
AML2.obs.loc[AML2.obs['orig.ident'].str.startswith("P"),'sam']=AML2.obs['orig.ident']
AML2.obs['sam']=AML2.obs['sam'].cat.remove_unused_categories()
AML2.obs['sam']=AML2.obs['sam'].cat.add_categories(AML2.obs['sample'].unique())
import re
r = re.compile('^P')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
b=AML2.obs['orig.ident'].unique()
b=b[vmatch(b)]
AML2.obs['sample']=AML2.obs['sample'].cat.add_categories(b)
AML2.obs.loc[AML2.obs['orig.ident'].str.startswith("N"),'sam']=AML2.obs['sample']


sc.pp.regress_out(AML2, ['total_counts', 'percent.mt','pct_counts_ribo'])
sc.pp.scale(AML2, max_value=10)
sc.tl.pca(AML2, svd_solver='arpack')
sc.pl.pca_variance_ratio(AML2, log=True, n_pcs = 50)

sce.pp.scanorama_integrate(AML2, key = 'sam',knn=40,sigma=20)
sc.pp.neighbors(AML2,n_pcs = 40, n_neighbors = 20)

sc.tl.leiden(AML2, resolution = 0.89, key_added = "leiden_089")
sc.tl.umap(AML2)
sc.pl.umap(AML2, color=['leiden_089'],save='umap_leiden_089_0502.png',legend_loc='on data')
sc.pl.umap(AML2, color='state', save='state_umap.png')
sc.pl.umap(AML2, color='sample', save='sample_umap.png')
marker_genes = ['CD79A','MS4A1',"CD19",'FCER1G','FCGR3A','CD14','KLF1','SNCA','HBD',"MPO","ELANE","HOPX","CD163","GATA1","GNLY","NKG7","CD27",'MZB1','PPBP','PF4',"GP9","CD34","CD38","IGLL1","MME","CD3D","CD3E","CD3G","CD1C","IL3RA","IRF8",'FCER1A','MSI2',"CDK6","CTSG"]
sc.pl.dotplot(AML2, marker_genes, groupby='leiden_089');

sc.tl.rank_genes_groups(AML2, 'leiden_089', method='t-test', key_added = "t-test089")
pd.DataFrame(AML2.uns['t-test05']['names']).head(5)
result = AML2.uns['t-test089']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).to_csv('/data/AMLscRNA/mergeblood/20220501/marker_089dedouble.csv', index=False,header=True)

cluster2annotation = {
'0':'undifferentiated1',
'1':'undifferentiated2',
'2':'GMP',
'3':'T',
'4':'undifferentiated3',
'5':'undifferentiated4',
'6':'HSPC',
'7':'undifferentiated5',
'8':'undifferentiated6',
'9':'CD14+mono',
'10':'MEP',
'11':'undifferentiated7',
'12':'undifferentiated2',
'13':'undifferentiated8',
'14':'Early_ery',
'15':'T',
'16':'undifferentiated9',
'17':'NK',
'18':'cDC',
'19':'Early_ery',
'20':'undifferentiated9',
'21':'Ery',
'22':'CD14+mono',
'23':'pDC',
'24':'B',
'25':'CD16+mono',
'26':'M2_macro',
'27':'undifferentiated5',
'28':'undifferentiated8',
'29':'undifferentiated1',
'30':'Pre-B',
'31':'Plasma',
'32':'Platelets',
'33':'Stromal cell',
'34':'undifferentiated1',
'35':'undifferentiated1',
}

AML2.obs['celltype'] = AML2.obs['leiden_089'].map(cluster2annotation).astype('category')
import seaborn as sns
sns.set_style('whitegrid')
sns.set(rc={"figure.figsize":(8, 8)})
sc.pl.umap(AML2, color='celltype', legend_loc='on data',frameon=False, legend_fontsize=8,size=2, legend_fontoutline=0,save='umap_leiden_089_celltype.pdf')
import re
r = re.compile('P[1-6]b')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
AML2.obs['type']=AML2.obs['type'].cat.add_categories(["RESb","RESa","UNb","UNa","HD"])
AML2.obs.loc[vmatch(AML2.obs['orig.ident']),'type']="RESb"
r = re.compile('P[1-6]a')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
AML2.obs.loc[vmatch(AML2.obs['orig.ident']),'type']="RESa"
r = re.compile('P[7-9]b|P10b')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
AML2.obs.loc[vmatch(AML2.obs['orig.ident']),'type']="UNb"
r = re.compile('P[7-9]a|P10a')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
AML2.obs.loc[vmatch(AML2.obs['orig.ident']),'type']="UNa"
AML2.obs.loc[AML2.obs['batch'] == 'health','type']="HD"

AML2.write_h5ad('AML2_annotate0501.h5ad')
marker_genes = ['CD79A','MS4A1',"CD19",'FCER1G','FCGR3A','CD14','KLF1','SNCA','HBD',"MPO","ELANE","HOPX","CD34","CD163","GATA1","GNLY","NKG7","KLRF1","CD27",'MZB1','PPBP','PF4',"GP9","CD38","IGLL1","MME",'MGP','LEPR',"CD3D","CD3E","CD3G","CD1C","IL3RA","IRF8",'FCER1A','MSI2',"CDK6","CTSG"]
sns.set(rc={"figure.figsize":(12, 8)}) 
sns.set_style('ticks')
sc.pl.dotplot(AML2, marker_genes, groupby='celltype',save='leiden_089_celltypedot.pdf');
pd.crosstab(columns=AML2.obs['celltype'],index =AML2.obs['orig.ident']).to_csv('sample_celltype0429.csv', index=True,header=True,sep="\t")
my_palette=dict(zip(AML2.obs['celltype'].sort_values().unique(),['#023fa5', '#7d87b9', '#bec1d4', '#d6bcc0', '#bb7784', '#8e063b','#4a6fe3', '#8595e1', '#b5bbe3', '#e6afb9', '#e07b91', '#d33f6a','#11c638', '#8dd593', '#c6dec7', '#ead3c6', '#f0b98d', '#ef9708','#0fcfc0', '#9cded6', '#d5eae7', '#f3e1eb', '#f6c4e1', '#f79cd4','#7f7f7f', '#c7c7c7', '#1ce6ff']))
color_palette={'B': '#023fa5', 'CD14+mono': '#7d87b9', 'CD16+mono': '#bec1d4', 'Early_ery': '#d6bcc0', 'Ery': '#bb7784', 'GMP': '#8e063b', 'HSPC': '#4a6fe3', 'M2_macro': '#8595e1', 'MEP': '#b5bbe3', 'NK': '#e6afb9', 'Plasma': '#e07b91', 'Platelets': '#d33f6a', 'Pre-B': '#11c638', 'Stroma': '#8dd593', 'T': '#c6dec7', 'cDC': '#ead3c6', 'pDC': '#f0b98d'}
heal=AML2[~AML2.obs['celltype'].str.startswith("undiff")]
import scanyuan as scy
sns.set(rc={"figure.figsize":(20, 8)}) 
sns.set(style="ticks")
scy.stacked_violin_t(heal, marker_genes, size=0.5,groupby='celltype',save='celltype089_violin.pdf',palette=my_palette)

dfa = sc.get.obs_df(P9, ["BCL2L1", 'state', 'CellTypeN'], use_raw=False)
df1 = dfa.melt(id_vars=["CellTypeN",'state'], value_vars="BCL2L1")
import seaborn as sns
sns.set(rc={"figure.figsize":(16, 16)}) 
sns.set(style="ticks")
g=sns.catplot(x = "state", y = "value", hue = "CellTypeN",sharey=False,linewidth=0.4, kind = 'violin',order=order, col="CellTypeN", data = df1,palette = keyy,col_order=sort1,height=3,aspect=1,col_wrap=6,sharex=True,inner="quartile",hue_order=groupby_order,dodge=False)
g.set_xticklabels(rotation=45, fontsize = 12)
g.set_yticklabels(fontsize = 12)
for (row_key),ax in g.axes_dict.items():
    ax.set(ylabel=None)
    ax.set_title(f"{row_key}", horizontalalignment='center',size=16)
    ax.tick_params(which='major',direction='in',length=2,width=0.5,labelsize=12)

g.tight_layout(pad=2, w_pad=2, h_pad=2) 
g.set(xlabel=None)
g.fig.subplots_adjust(hspace=0.3) 
g.despine(top=False,right=False,left=False,bottom=False)
import matplotlib.pyplot as plt
plt.savefig("celltype_BCL2L1_allgene.pdf")

#AML cell
a=AML2.obs['celltype'].str.startswith("undiff")&AML2.obs['sample'].str.startswith("P")
a = a[a.index.str.startswith("P")]
undiff = AMLy[a,:]
undiff.write_h5ad('undiff0505.h5ad')


#correlation matrix
   
for clust in AML2.obs['CellTypeN'].cat.categories: 
    res.loc[clust] = AML2[AML2.obs['CellTypeN'].isin([clust]),:].X.mean(0)
X_sim=res.transpose()
X_sim.to_csv('/data/AMLscRNA/mergeblood/20220501/for_correlation.csv', index=True,header=True,sep="\t")
list(X_sim)
X_sim= pd.DataFrame(X_sim, columns=['B', 'CD14+mono', 'CD16+mono', 'Early_ery', 'Ery', 'GMP', 'GMP_like', 'HSPC', 'LSPC_Cycle', 'LSPC_Primed', 'LSPC_Quiescent', 'MEP', 'Mono_like', 'NK', 'Plasma', 'Platelets', 'Pre-B', 'ProMono_like', 'Stroma', 'T', 'cDC', 'cDC_like', 'pDC'])
my_r = np.corrcoef(X_sim)
'''
cell.prop2<-read.csv("/data/AMLscRNA/mergeblood/20220501/correlation_matrix1",header = TRUE,sep="\t",row.names="type")
col_fun = colorRamp2(c(0.42,0.7,0.99), c("#1F77B4", "white", "#D52728"))
Heatmap(cell.prop2, name=" ", col=col_fun,border=TRUE,rect_gp=gpar(col="white",lwd=2))
'''

AML2.raw=sc.read_h5ad("/data/AMLscRNA/mergeblood/AML10pC_blood5n.h5ad")
AML10=AML2[AML2.obs['batch'] == 'patient',:]
AML10.obs.index = [i.split('-')[0]+"-1" for i in AML10.obs.index]
celltype=pd.read_csv("/data/AMLscRNA/mergeblood/20220501/undiff/scArches_PVG_DAV_CellType_0501.csv", delimiter=',', index_col=0)
celltype.index = [i.split('-')[0]+"-1" for i in celltype.index]
c=pd.merge(AML10.obs,celltype,how='left', left_index=True, right_index=True)
c=c.drop(['Cohort','batch_y','state_y'],axis=1)
c.rename(columns = {'predictions':'CellTypeN'}, inplace = True)
c.loc[c['CellTypeN'].isnull(),'CellTypeN']=c['celltype']
AML10.obs=c
AML10.obs.rename(columns = {'state_x':'state','batch_x':'batch'}, inplace = True)
AML10.write_h5ad('/data/AMLscRNA/mergeblood/20220501/AML10_0502.h5ad')

AML10=sc.read_h5ad("/data/AMLscRNA/mergeblood/20220501/AML10_0502.h5ad")
comp=AML10[AML10.obs['celltype'].str.startswith("undiff"),:]
marker_genes=['BCL2','MCL1','KIT','ITGAM','CD68','FCGR1A']
marker_genes=["BRD2","BRD4","CSF1R","NFKB1","PIK3CD","SYK","VEGFA"]
keyy={'B': '#023fa5', 'CD14+mono': '#7d87b9', 'CD16+mono': '#bec1d4', 'Early_ery': '#d6bcc0', 'Ery': '#bb7784', 'GMP': '#8e063b', 'HSPC': '#4a6fe3', 'M2_Macro': '#8595e1', 'MEP': '#b5bbe3', 'NK': '#e6afb9',  'Plasma': '#d33f6a', 'Platelets': '#11c638', 'Pre-B': '#8dd593', 'T': '#c6dec7', 'cDC': '#ead3c6', 'pDC': '#f0b98d', 'Stroma': '#ef9708','LSPC_Cycle':'#FF7F0E', 'LSPC_Primed':'#279E68', 'LSPC_Quiescent':'#D62728','Mono_like':'#AA40FC', 'ProMono_like':'#8C564B', 'cDC_like':'#E377C2', 'GMP_like':'#1F77B4'}
groupby_order = ['HSPC',  'GMP','MEP','Pre-B','B','Plasma', 'T', 'NK',  'CD14+mono', 'CD16+mono','M2_Macro','pDC','cDC','Early_ery','Ery','Platelets','Stroma',  'LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
groupby_order = ['LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
C = {val:ix for ix,val in enumerate(groupby_order)}	
dfa = sc.get.obs_df(comp, ['CellTypeN',*marker_genes], use_raw=True)
df1 = dfa.melt(id_vars=["CellTypeN"], value_vars=marker_genes)
sort1=sorted(df1['CellTypeN'].unique().tolist(), key=C.get)	
sns.set(style="ticks")
g=sns.catplot(x = "CellTypeN", y = "value", hue = "variable", kind = 'violin', col="variable", cut=0,data = df1,height=1,aspect=20,col_wrap=1,inner=None,sharex=True,order=sort1,hue_order=marker_genes,dodge=False)
g.set_xticklabels(rotation=45)
for (row_key),ax in g.axes_dict.items():
    ax.set_ylabel(f"{row_key}", rotation=0, horizontalalignment='right', fontsize = 15)

ax.set_xticklabels(CellTypeN, fontsize=15) 
g.tight_layout(pad=1, w_pad=1, h_pad=1) 
g.set(title=None)
g.fig.subplots_adjust(hspace=0) 
g.despine(top=False,right=False,left=False,bottom=False)
plt.savefig("celltype_gene.pdf")
plt.savefig("celltype_gene_tumor.pdf")

marker_genes=['BCL2','BCL2L1','BCL2L2','MCL1','BCL2L10','BAX','BAK1','BID','BCL2L11','PMAIP1','BBC3','BMF','BAD','BIK','HRK']
keyy={'B': '#023fa5', 'CD14+mono': '#7d87b9', 'CD16+mono': '#bec1d4', 'Early_ery': '#d6bcc0', 'Ery': '#bb7784', 'GMP': '#8e063b', 'HSPC': '#4a6fe3', 'M2_Macro': '#8595e1', 'MEP': '#b5bbe3', 'NK': '#e6afb9',  'Plasma': '#d33f6a', 'Platelets': '#11c638', 'Pre-B': '#8dd593', 'T': '#c6dec7', 'cDC': '#ead3c6', 'pDC': '#f0b98d', 'Stroma': '#ef9708','LSPC_Cycle':'#FF7F0E', 'LSPC_Primed':'#279E68', 'LSPC_Quiescent':'#D62728','Mono_like':'#AA40FC', 'ProMono_like':'#8C564B', 'cDC_like':'#E377C2', 'GMP_like':'#1F77B4'}
groupby_order = ['HSPC',  'GMP','MEP','Pre-B','B','Plasma', 'T', 'NK',  'CD14+mono', 'CD16+mono','M2_Macro','pDC','cDC','Early_ery','Ery','Platelets','Stroma',  'LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
C = {val:ix for ix,val in enumerate(groupby_order)}	
dfa = sc.get.obs_df(AML10, ['CellTypeN',*marker_genes], use_raw=True)
df1 = dfa.melt(id_vars=["CellTypeN"], value_vars=marker_genes)
sort1=sorted(df1['CellTypeN'].unique().tolist(), key=C.get)	
sns.set(style="ticks")
g=sns.catplot(x = "CellTypeN", y = "value", hue = "variable", kind = 'violin', col="variable", cut=0,data = df1,height=1,aspect=20,col_wrap=1,inner=None,sharex=True,order=sort1,hue_order=marker_genes,dodge=False)
g.set_xticklabels(rotation=45)
for (row_key),ax in g.axes_dict.items():
    ax.set_ylabel(f"{row_key}", rotation=0, horizontalalignment='right', fontsize = 15)

g.tight_layout(pad=1, w_pad=1, h_pad=1) 
g.set(title=None)
g.fig.subplots_adjust(hspace=0)
g.despine(top=False,right=False,left=False,bottom=False)
plt.savefig("celltype_BCL2gene.pdf")

comp=AML10[(AML10.obs['celltype'].str.startswith("undiff"))&(AML10.obs['state'] == 'before'),:]
comp=AML10[(AML10.obs['celltype'].str.startswith("undiff")|(AML10.obs['celltype'] == 'Stromal cell'))&(AML10.obs['state'] == 'before'),:]
comp=heal[(heal.obs['celltype']=="HSPC")&((heal.obs['state'] == 'before')|(heal.obs['type'] == 'HD')),:]
marker_genes=['KIT','CD14','CD34','MAFB']
y="KIT"
y="CD34"
y="THY1"
y="CD14"
y="MAFB"
y="CD68"
y="ITGAM"
y="SOD2"
y="GRK2"
y="BHLHE40"
y="SMAD1"
y="JUNB"
y="CXCR4"
y="CXCL12"
y="KLF8"
dfa = sc.get.obs_df(comp, ['CellTypeN',y,'type'], use_raw=True)
df1 = dfa.melt(id_vars=["CellTypeN",'type'], value_vars=[y])
groupby_order = ['LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
order = ['RESb','UNb']
C = {val:ix for ix,val in enumerate(groupby_order)}	
sort1=sorted(df1['CellTypeN'].unique().tolist(), key=C.get)	
sns.set(style="ticks")
groupby_order = ['LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like','Stroma']
g=sns.catplot(x = "type", y = "value", hue = "CellTypeN", kind = 'violin', col="CellTypeN", cut=0, data = df1,height=3,aspect=0.5,col_wrap=8,inner=None,sharey=True,order=order,col_order=groupby_order,dodge=False)
g.set_xticklabels(rotation=45)
g.set_ylabels(y, rotation=90)
for (row_key),ax in g.axes_dict.items():
    ax.set_title(f"{row_key}", rotation=0, fontsize = 12)
    ax.set_xlabel(xlabel = None)
       
ax.set_xticklabels(order, fontsize=12) 
g.tight_layout(pad=1, w_pad=1, h_pad=1) 
g.fig.subplots_adjust(wspace=0)
g.despine(top=False,right=False,left=False,bottom=False)
plt.savefig(y+"_gene.pdf")

comp=AML10[(AML10.obs['celltype'].str.startswith("undiff"))&(AML10.obs['state'] == 'before'),:]
marker_genes=['BCL2','BCL2L1','BCL2L2','MCL1','BCL2L10','BAX','BAK1','BID','BCL2L11','PMAIP1','BBC3','BMF','BAD','BIK','HRK']
marker_genes=["BRD2","BRD4","CSF1R","NFKB1","PIK3CD","SYK","VEGFA"]
y="BCL2"
y="BCL2L1"
y="BCL2L2"
y="BCL2A1"
y="MCL1"
y="BCL2L10"
y="BAX"
y="BAK1"
y="BID"
y="BCL2L11"
y="PMAIP1"
y="BBC3"
y="BMF"
y="BAD"
y="BIK"
y="HRK"
y="AXL"
y="CLPB"
y="CDK9"
y="CXCR4"
y="CD44"

comp=AML10[(AML10.obs['celltype'].str.startswith("undiff"))&(AML10.obs['state']=="before"),:]
y="BRD2"
y="BRD4"
y="CSF1R"
y="NFKB1"
y="PIK3CD"
y="SYK"
y="VEGFA"
y="MCL1"
y="PTPN6"
y="ATF4"
dfa = sc.get.obs_df(comp, ['CellTypeN',y,'type'], use_raw=True)
df1 = dfa.melt(id_vars=["CellTypeN",'type'], value_vars=[y])
groupby_order = ['HSPC',  'GMP','MEP','Pre-B','B','Plasma', 'T', 'NK',  'CD14+mono', 'CD16+mono','M2_Macro','pDC','cDC','Early_ery','Ery','Platelets','Stroma',  'LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
order=['RESb','RESa','UNb','UNa']
groupby_order = ['LSPC_Primed','LSPC_Quiescent','LSPC_Cycle', 'GMP_like', 'ProMono_like', 'cDC_like', 'Mono_like']
order = ['RESb','UNb']
C = {val:ix for ix,val in enumerate(groupby_order)}	
sort1=sorted(df1['CellTypeN'].unique().tolist(), key=C.get)	
import seaborn as sns
sns.set(style="ticks")
g=sns.catplot(x = "type", y = "value", hue = "type", kind = 'violin', col="CellTypeN", cut=0, data = df1,height=2.5,aspect=1,col_wrap=7,inner=None,sharey=True,order=order,col_order=sort1,dodge=False)
g.set_ylabels(y, rotation=90)
for (row_key),ax in g.axes_dict.items():
    ax.set_title(f"{row_key}", rotation=0, fontsize = 12)
    ax.set_xlabel(xlabel = None)
    ax.set_xticklabels(labels=order,rotation=45)
       
ax.set_xticklabels(order, fontsize=12) 
g.tight_layout(pad=1, w_pad=1, h_pad=1) 
g.fig.subplots_adjust(wspace=0) 
g.fig.subplots_adjust(hspace=0) 
g.despine(top=False,right=False,left=False,bottom=False)
plt.savefig(y+"_geneT.pdf")

comp=heal[(heal.obs['celltype']=="HSPC")&((heal.obs['state'] == 'before')|(heal.obs['type'] == 'HD')),:]
marker_genes=['BCL2','MCL1','KIT','CD14','CD34','MAFB']
comp=AML10[(AML10.obs['celltype']=="HSPC"),:]
comp=AML10[AML10.obs['celltype'].str.startswith("undiff"),:]
marker_genes=['MDK','MLLT3','PFDN5','RACK1','STK4']
marker_genes=['S100A9','S100A8','NFKBIA','HSPA1B','IL1B','CLU','RPS27A','HSPA1A','RPS3','UBC','RIPK2']
y="KIT"
y="CD14"
y="CD34"
y="MAFB"
dfa = sc.get.obs_df(comp, [*marker_genes,'type'], use_raw=True)
df1 = dfa.melt(id_vars=['type'], value_vars=marker_genes)
order=['RESb','RESa','UNb','UNa']
C = {val:ix for ix,val in enumerate(groupby_order)}	
sort1=sorted(df1['CellTypeN'].unique().tolist(), key=C.get)	
sns.set(style="ticks")
g=sns.catplot(x = "type", y = "value", hue = "variable", kind = 'violin', col="variable", cut=0, data = df1,height=3,aspect=0.5,col_wrap=6,inner=None,sharey=True,sharex=True,order=order,dodge=False)
g.set_xticklabels(rotation=45)
g.set_ylabels("unresponder HSPC", rotation=90)
for (row_key),ax in g.axes_dict.items():
    ax.set_title(f"{row_key}", rotation=0, fontsize = 12)
    ax.set_xlabel(xlabel = None)
       
ax.set_xticklabels(order, fontsize=12) 
g.tight_layout(pad=1, w_pad=1, h_pad=1)  
g.fig.subplots_adjust(wspace=0) 
g.fig.subplots_adjust(hspace=0) 
g.despine(top=False,right=False,left=False,bottom=False)
plt.savefig("HSPC_wntgene.pdf")
plt.savefig("HSPC_NFKBgene.pdf")

#diff gene and GSEA
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(anndata)
res = pd.DataFrame(columns=AML2.obs['CellTypeN'], index=AML2.obs.index)   
AML10.obs.to_csv('/data/AMLscRNA/mergeblood/20220501/AML10_obs.csv',index=True,header=True)
AML10.obs<-read.csv('/data/AMLscRNA/mergeblood/20220501/AML10_obs.csv',row.names = "X")
c<-merge(AML10.obs,AML@meta.data,by=0)
c=c[,-c(27:34)]
colnames(c)<-gsub(".x","",colnames(c))
rownames(c)=c[,1] 
c=c[,-1]
AML@meta.data<-c
AML$celltype.condi <- paste(AML$type,AML$CellTypeN,sep = "_")
humandata<- AML
Idents(AML) <- "celltype.condi"
humandata <- NormalizeData(humandata, normalization.method = "LogNormalize", scale.factor = 10000)

humandata@meta.data<-AML@meta.data
humandata@meta.data$act<-gsub('.{1}$', '',humandata@meta.data$type)

LSPC_Quiescent1=subset(humandata,subset=(CellTypeN =="LSPC_Quiescent"&&act =="RES"))
Idents(LSPC_Quiescent1) <- "celltype.condi"
LSPC_Quiescent_RES<-FindMarkers(LSPC_Quiescent1, ident.1 = "RESa_LSPC_Quiescent", ident.2 = "RESb_LSPC_Quiescent",test.use = "MAST")

LSPC_Quiescent_RES<-FindMarkers(humandata, ident.2 = "RESb_LSPC_Quiescent", ident.1 = "RESa_LSPC_Quiescent",test.use = "MAST")
LSPC_Quiescent_RES_NB<-FindMarkers(humandata, ident.2 = "RESb_LSPC_Quiescent", ident.1 = "RESa_LSPC_Quiescent",test.use = "negbinom",slot="counts")
write.table(LSPC_Quiescent_RES_NB,"/data/AMLscRNA/mergeblood/20220501/LSPC_Quiescent/LSPC_Quiescent_RES_NB.bed",quote = F,sep="\t")
LSPC_Quiescent=subset(AML,subset=(CellTypeN =="LSPC_Quiescent"))
LSPC_Quiescent1=subset(humandata,subset=(CellTypeN =="LSPC_Quiescent"))
LSPC_Quiescent1$celltype.condi <- paste(LSPC_Quiescent1$type,LSPC_Quiescent1$CellTypeN,sep = "_")
Idents(LSPC_Quiescent1) <- "celltype.condi"
VlnPlot(object = LSPC_Quiescent, features = c('APOC2','ANXA2R'),  pt.size = 0, combine = FALSE, log = TRUE)
VlnPlot(object = LSPC_Quiescent1, features = c('BIK'),  pt.size = 0.2, combine = FALSE, log = FALSE)
RidgePlot(object = LSPC_Quiescent1, features = c('AZU1'))
all.markers = diff.mast %>% select(gene, everything()) %>% subset(p_val<0.05)
write.table(LSPC_Quiescent_RES,"/data/AMLscRNA/mergeblood/20220501/LSPC_Quiescent/LSPC_Quiescent_RES2.bed",quote = F,sep="\t")
saveRDS(humandata,"/data/AMLscRNA/mergeblood/20220501/Rdata_forDEG.rds.gz",compress="gz")

rnk = pd.read_csv("/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC-HSPC.rnk", header=None, sep="\t")
pre_res = gseapy.prerank(rnk=rnk, gene_sets='GO_Biological_Process_2021',processes=4,permutation_num=100, outdir='/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC-HSPC_goBP', format='png', seed=6)
rnk = pd.read_csv("/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC_Primed-HSPC.rnk", header=None, sep="\t")
pre_res = gseapy.prerank(rnk=rnk, gene_sets='GO_Biological_Process_2021',processes=4,permutation_num=100, outdir='/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC_Primed-HSPC_goBP1', format='png', seed=6,graph_num=60)
rnk = pd.read_csv("/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC_Cycle-HSPC.rnk", header=None, sep="\t")
pre_res = gseapy.prerank(rnk=rnk, gene_sets='GO_Biological_Process_2021',processes=4,permutation_num=100, outdir='/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC_Cycle-HSPC_goBP', format='png', seed=6)
rnk = pd.read_csv("/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC_Quiescent-HSPC.rnk", header=None, sep="\t")
pre_res = gseapy.prerank(rnk=rnk, gene_sets='GO_Biological_Process_2021',processes=4,permutation_num=100, outdir='/data/AMLscRNA/mergeblood/20220501/MAST/compare/GSEA/LSPC_Quiescent-HSPC_goBP', format='png', seed=6)
gseaplot(rank_metric=pre_res.ranking, term=terms[0], **pre_res.results[terms[0]])


rnk = pd.read_csv("/data/AMLscRNA/mergeblood/20220501/MAST/GSEA/HSPC_UN_diff.rnk", header=None, sep="\t")
pre_res = gseapy.prerank(rnk=rnk, gene_sets='GO_Biological_Process_2021',processes=4,permutation_num=100, outdir='/data/AMLscRNA/mergeblood/20220501/MAST/GSEA/HSPC_UN1', format='png', seed=6,graph_num=125)


#cytotrace
'''
counts <-  as.matrix(GetAssayData(humandata, slot = "data"))
write.csv(counts,"/data/AMLscRNA/mergeblood/20220501/cytotrace/matrix_normalize.csv",quote = F)
cat /data/AMLscRNA/mergeblood/20220501/AML10_obs.csv|awk -F ',' '{print$1"\t"$21}'|sed '1d'>/data/AMLscRNA/mergeblood/20220501/cytotrace/batch
cat /data/AMLscRNA/mergeblood/20220501/AML10_obs.csv|awk -F ',' '{print$1"\t"$NF}'|sed '1d'>/data/AMLscRNA/mergeblood/20220501/cytotrace/phenotype
'''


AML2.obs['after'] = pd.Categorical((AML2.obs['state'] == "after"))
AML2.uns['after_colors'] = ["#CECCCD", "#FF7F0E"]
sc.pl.umap(AML2, color='after', save='after_state.png')

AML2.obs['HD'] = pd.Categorical((AML2.obs['type'] == "HD"))
AML2.uns['HD_colors'] = ["#CECCCD", "#279E68"]
sc.pl.umap(AML2, color='HD', save='HD_state.png')

AML2.obs['before'] = pd.Categorical((AML2.obs['state'] == "before"))
AML2.uns['before_colors'] = ["#CECCCD", "#D62728"]
sc.pl.umap(AML2, color='before', save='before_state.png')

AML2.obs['RESb'] = pd.Categorical((AML2.obs['type'] == "RESb"))
AML2.uns['RESb_colors'] = ["#CECCCD", "#17BECF"]
sc.pl.umap(AML2, color='RESb', save='RESb_state.png')

AML2.obs['RESa'] = pd.Categorical((AML2.obs['type'] == "RESa"))
AML2.uns['RESa_colors'] = ["#CECCCD", "#98DF8A"]
sc.pl.umap(AML2, color='RESa', save='RESa_state.png')

AML2.obs['UNb'] = pd.Categorical((AML2.obs['type'] == "UNb"))
AML2.uns['UNb_colors'] = ["#CECCCD", "#FF9896"]
sc.pl.umap(AML2, color='UNb', save='UNb_state.png')

AML2.obs['UNa'] = pd.Categorical((AML2.obs['type'] == "UNa"))
AML2.uns['UNa_colors'] = ["#CECCCD", "#FFBB78"]
sc.pl.umap(AML2, color='UNa', save='UNa_state.png')







AML=sc.read_h5ad('/data/AMLscRNA/mergeblood/20220501/AML10_with_raw.h5ad')
raw=AML.raw.to_adata()
raw=raw[(raw.obs['CellTypeN'] != 'Ery')&(raw.obs['CellTypeN'] != 'Early_ery'),:]
RESb=raw[raw.obs['type'] == 'RESb',:]
RESa=raw[raw.obs['type'] == 'RESa',:]
UNb=raw[raw.obs['type'] == 'UNb',:]
UNa=raw[raw.obs['type'] == 'UNa',:]
RESb.write_h5ad('/data/AMLscRNA/mergeblood/20220501/cellchat/RESb.h5ad')
RESa.write_h5ad('/data/AMLscRNA/mergeblood/20220501/cellchat/RESa.h5ad')
UNb.write_h5ad('/data/AMLscRNA/mergeblood/20220501/cellchat/UNb.h5ad')
UNa.write_h5ad('/data/AMLscRNA/mergeblood/20220501/cellchat/UNa.h5ad')
 