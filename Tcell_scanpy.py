import numpy as np
import pandas as pd
import scanpy as sc
import sc_toolbox.api as sct
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 #保证输出pdf文字可以编辑
import scanpy.external as sce
import matplotlib.pyplot as plt
import os
os.chdir('/data/AMLscRNA/mergeblood/20220501/Tcell')
AML = sc.read_h5ad('/data/AMLscRNA/mergeblood/AML10pC_blood5n.h5ad') 
AML2 = sc.read_h5ad('/data/AMLscRNA/mergeblood/20220501/AML2_annotate0501.h5ad') #一共132025
Tcell = AML[(AML2.obs['celltype'] == 'T'),:]
Tcell.write_h5ad("/data/AMLscRNA/mergeblood/20220501/Tcell/Tcell0501.h5ad")
Tcell=sc.read_h5ad("/data/AMLscRNA/mergeblood/20220501/Tcell/Tcell0501.h5ad")
# ribosomal genes
Tcell.var['ribo'] = Tcell.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
Tcell.var['hb'] = Tcell.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(Tcell, qc_vars=['ribo','hb'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(Tcell, ['n_genes_by_counts', 'total_counts', 'percent.mt','pct_counts_ribo', 'pct_counts_hb'],jitter=0.4, groupby = 'orig.ident', rotation= 45,size=0.1,save='Tcell_QC_dedoublet.png')
sc.pl.highest_expr_genes(Tcell, n_top=20,save='Tcell_Top20_dedoublet.png')
sc.pl.violin(Tcell, ['percent.mt'],jitter=0.4, groupby = 'leiden_1', rotation= 45,size=0.1,save='Tcell_QC_dedoublet.png')
sc.pl.violin(T2, ['percent.mt'],jitter=0.4, groupby = 'leiden_2', rotation= 45,size=0.1,save='Tcell_QC_dedoublet.png')

mito_genes = Tcell.var_names.str.startswith('MT-')
Tcell.obs['percent_mt2'] = np.sum(
    Tcell[:, mito_genes].X, axis=1).A1 / np.sum(Tcell.X, axis=1).A1
Tcell.obs['n_counts'] = Tcell.X.sum(axis=1).A1

sc.pl.scatter(Tcell, x='total_counts', y='percent.mt', color="orig.ident",save='Tcell_count_MT_scatter_dedoublet.png')

cell_cycle_genes = [x.strip() for x in open('/data/AMLscRNA/cell_score')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in Tcell.var_names]


sc.pp.normalize_total(Tcell, target_sum=1e4)
sc.pp.log1p(Tcell)
Tcell.raw = Tcell
sc.tl.score_genes_cell_cycle(Tcell, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(Tcell, ['S_score', 'G2M_score'],jitter=0.4, groupby = 'orig.ident', rotation=45,size=0.1,save='Tcell_cellcyle0407.png')


Tcell.obs.loc[Tcell.obs['orig.ident'].str.startswith("P"),'sam']=Tcell.obs['orig.ident']
Tcell.obs['sam']=Tcell.obs['sam'].cat.remove_unused_categories() 
Tcell.obs['sam']=Tcell.obs['sam'].cat.add_categories(Tcell.obs['sample'].unique())
import re
r = re.compile('^P')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
b=Tcell.obs['orig.ident'].unique()
b=b[vmatch(b)]
Tcell.obs['sample']=Tcell.obs['sample'].cat.add_categories(b) 
Tcell.obs.loc[Tcell.obs['orig.ident'].str.startswith("N"),'sam']=Tcell.obs['sample']
Tcell2 = Tcell.raw.to_adata() 
T=Tcell2
sc.pp.highly_variable_genes(Tcell2, min_mean=0.0125, max_mean=3, min_disp=0.5,subset=True)

sc.pp.regress_out(Tcell2, ['total_counts', 'percent.mt','pct_counts_ribo'])
sc.pp.scale(Tcell2, zero_center=False)
sc.tl.pca(Tcell2, svd_solver='arpack')
sc.pl.pca_variance_ratio(Tcell2, log=True, n_pcs = 50)#30

pd.crosstab(columns=Tcell2.obs['state'], index=Tcell2.obs['orig.ident'])

sce.pp.harmony_integrate(Tcell2, key = 'sam',knn=40,sigma=20)
sc.pp.neighbors(Tcell2, n_pcs = 40, n_neighbors = 20, use_rep="X_pca_harmony")

sc.tl.leiden(Tcell2, resolution = 2, key_added = "leiden_2")
sc.tl.umap(Tcell2)
sc.pl.umap(Tcell2, color=['leiden_2'],legend_loc='on data')
marker_genes = ['CD4','CD8A',"CCR7","SELL","LEF1","TCF7","CD69","GZMH","GNLY","NKG7","RORC","KLRB1","IL2RA","MKI67","PCNA","HAVCR2","PRDM1","TIGIT","TOX","IFIT3","ISG15","ISG20","MX1","IFNAR1","PRF1","PDCD1","RUNX3","KLRG1","ITGAE","B3GAT1","IL7R","IL4","CCR6","CCR10","IL2RA","GZMA","GZMB","GZMK","PRF1","IFNG","TRGC2","TRDC","FOXP3","CTLA4","IL15","LAG3","PDCD1"]
sc.pl.dotplot(Tcell2, marker_genes, groupby='leiden_2');

pd.crosstab(columns=Tcell2.obs['state'], index=T2.obs['leiden_2'])

Tcell.obs['RK']	='Keep'
Tcell.obs['ID'] = Tcell.obs.index.values.tolist()
remo=Tcell2[(Tcell2.obs['leiden_2'] == '8')|(Tcell2.obs['leiden_2'] == '15')|(Tcell2.obs['leiden_2'] == '16')|(Tcell2.obs['leiden_2'] == '17')|(Tcell2.obs['leiden_2'] == '22')|(Tcell2.obs['leiden_2'] == '28')|(Tcell2.obs['leiden_2'] == '29')|(Tcell2.obs['leiden_2'] == '32'),:].obs.index.values.tolist()#PPBP
anno = { i : j for i, j in zip(Tcell.obs['ID'],Tcell.obs['RK']) }
Tanno1 = { i : "remo" for i in remo }
anno.update(Tanno1)
Tcell.obs['RK'] = (Tcell.obs['ID'].map(anno).astype("category"))
pd.crosstab(columns=Tcell2.obs['leiden_2'], index=Tcell2.obs['batch'])
T3=Tcell[Tcell.obs['RK']=="Keep",:] 
sc.pp.normalize_total(T3, target_sum=1e4)
sc.pp.log1p(T3)
T3.raw = T3
sc.tl.score_genes_cell_cycle(T3, s_genes=s_genes, g2m_genes=g2m_genes)

T2 = T3.raw.to_adata()
sc.pp.highly_variable_genes(T2, n_top_genes=2500,batch_key="sam",subset=True)
sc.pp.regress_out(T2, ['total_counts', 'percent.mt','pct_counts_ribo'])
sc.pp.scale(T2,zero_center=False)
sc.tl.pca(T2, svd_solver='arpack')
sc.pl.pca_variance_ratio(T2, log=True, n_pcs = 50,save='Tcell_pc_0407_dedoublet.png')#30
sce.pp.scanorama_integrate(T2, key = 'sam',knn=40,sigma=20)
sc.pp.neighbors(T2, n_pcs = 40, n_neighbors = 20)
sc.tl.leiden(T2, resolution = 1.8, key_added = "leiden_18")
sc.tl.umap(T2)
sc.pl.umap(T2, color=['leiden_18'],save='Tcell_umap_compare_1.8dedoublet.png')
marker_genes = ['CD4','CD8A','CD8B','CCR7','LEF1','SELL','IL7R','STAT1','FOXP3','CTLA4','MX1','IFI44L','GZMK','CD27','LTB','GNLY','PRF1','GZMB','TRGV9','SLC4A10','PF4','HBB','LYZ','S100A8','S100A9','PPBP','SOX4','ELANE','TRDC','TRGC1']
sc.pl.dotplot(T2, marker_genes, groupby='leiden_18');

sc.tl.rank_genes_groups(T2, 'leiden_18', method='t-test', key_added = "ttest")
sc.pl.rank_genes_groups_dotplot(T2, n_genes=5, key="ttest", groupby="leiden_18")



Tannotation = {
'0':'CD8_CTL',
'1':'GD',
'2':'CD8_EM',
'3':'CD8_Naive',
'4':'CD8_EM',
'5':'CD4_Naive',
'6':'CD4_LTB',
'7':'CD8_EM',
'8':'MAIT',
'9':'CD4_LTB',
'10':'CD8_Naive',
'11':'CD4_LTB',
'12':'CD8_Naive',
'13':'CD8_Naive',
'14':'CD8_CTL',
'15':'CD8_STAT1',
'16':'CD8_CTL',
'17':'Treg',
'18':'CyclingCD8',
'19':'CD8_EM',
'20':'GD_GZMK',
'21':'CyclingGD',
'22':'Th',
}
T2.obs['celltype'] = T2.obs['leiden_18'].map(Tannotation).astype('category')
sc.pl.umap(T2, color='celltype', legend_loc='on data',frameon=False, legend_fontsize=12, legend_fontoutline=0,save='leiden_18_anno0501_Tcell.pdf')
marker_genes = ['CD4','CD8A','CD8B','CCR7','LEF1','SELL','IL7R','LTB','MX1','IFI44L','IFIT3','GZMK','STAT1','GNLY','PRF1','GZMB','TRGV9','SLC4A10','TRDC','TRGC2','MKI67','ICOS','FOXP3','CTLA4']
sc.pl.dotplot(T2, marker_genes, groupby='celltype',save='Tcelltypedot_2.pdf');
pd.crosstab(columns=T2.obs['celltype'], index=T2.obs['sam']).to_csv('/data/AMLscRNA/mergeblood/20220501/Tcell/Tcell_celltype.csv', index=True,header=True,sep="\t")
T2.write_h5ad('T_annotate.h5ad')
my_palette=dict(zip(T2.obs['celltype'].sort_values().unique(),['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b','#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a','#ff9896']))

cytotoxic_gene=['KLRG1','IL7R','GNLY','SELL','TYROBP','TRDC','SAMD3','FCRL6','CCND3','BIN2','ARL4C','EMP3','SORL1','C1orf162','STK38','EFHD2','FAM65B','S1PR1','CD300A','SPON2','KLRF1','FGR','PLEK','TGFBR3','C1orf21','PLAC8','FCGR3A','S1PR5','FGFBP2','CX3CR1']

sc.tl.score_genes(T2, cytotoxic_gene, ctrl_size=len(cytotoxic_gene), gene_pool=None, n_bins=25, score_name='cytotoxic_score', random_state=0, copy=False, use_raw=None)
sc.pl.violin(T2, ['cytotoxic_score'],jitter=0.4, groupby = 'com', rotation= 45,size=0,save='cytotoxic_score_sam.pdf',order=groupby_order)

T2.obs.loc[T2.obs['batch'].str.startswith("h"),'com']=T2.obs['batch']
T2.obs['com']=T2.obs['com'].cat.add_categories(["RESb","RESa","UNb","UNa"]) 

import re
r = re.compile('P[1-6]b')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
T2.obs.loc[vmatch(T2.obs['orig.ident']),'com']="RESb"
r = re.compile('P[1-6]a')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
T2.obs.loc[vmatch(T2.obs['orig.ident']),'com']="RESa"
r = re.compile('P[7-9]b|P10b')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
T2.obs.loc[vmatch(T2.obs['orig.ident']),'com']="UNb"
r = re.compile('P[7-9]a|P10a')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
T2.obs.loc[vmatch(T2.obs['orig.ident']),'com']="UNa"
T2.obs.loc[T2.obs['batch']=="health",'com']="HD"
pd.DataFrame([T2.obs['orig.ident'],T2.obs['com']]).transpose().to_csv('/data/AMLscRNA/mergeblood/Tcell/check', index=False,header=True,sep="\t")
cytotoxic=T2[T2.obs['celltype'] == 'GD',:]
sc.tl.score_genes(cytotoxic, cytotoxic_gene, ctrl_size=len(cytotoxic_gene), gene_pool=None, n_bins=25, score_name='cytotoxic_score', random_state=0, copy=False, use_raw=None)
sc.pl.violin(cytotoxic, ['cytotoxic_score'],jitter=0.4, groupby = 'com', rotation= 45,size=0,save='cytotoxic_score_sam.pdf',order=groupby_order)
dsf=T2[T2.obs['celltype'] == 'CD8_exhausted',:]
sc.tl.score_genes(dsf, exhausted_gene, ctrl_size=len(exhausted_gene), gene_pool=None, n_bins=25, score_name='exhausted_score', random_state=0, copy=False, use_raw=None)

#GSVA
T=sc.read_h5ad("/data/AMLscRNA/mergeblood/20220429/Tcell/Tcell0429.h5ad")
Tcell=sc.read_loom("/data/AMLscRNA/mergeblood/Tcell0412.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
T2 = sc.read_h5ad('/data/AMLscRNA/mergeblood/20220501/Tcell/T2_annotate.h5ad') #一共132025
sc.pp.normalize_total(T, target_sum=1e4)
T.obs['ID'] = T.obs.index.values.tolist()
T.obs['RK']	='del'
anno = { i : j for i, j in zip(T.obs['ID'],T.obs['RK']) }
rem=T2.obs.index.values.tolist()
Tanno1 = { i : "Keep" for i in rem }
anno.update(Tanno1)
T.obs['RK'] = (T.obs['ID'].map(anno).astype("category"))
T3=Tcell[Tcell.obs['RK']=="Keep",:] 
from GSVA import gsva, gmt_to_dataframe
genesets_df = gmt_to_dataframe('/data/AMLscRNA/mergeblood/Tcell/temp/gs.gmt')
expression_df =pd.DataFrame(T3.to_df().transpose())
expression_df=expression_df.transpose()
expression_df= expression_df[(expression_df.T != 0).any()]
pathways_df = GSVA.gsva(expression_df,genesets_df,tempdir ="/data/AMLscRNA/mergeblood/20220429/Tcell/",ssgsea_norm=True,verbose=True)

GSVA=pd.read_csv("/data/AMLscRNA/mergeblood/20220501/Tcell/pathways1.csv", delimiter=',',index_col=0).transpose()
c=pd.merge(T2.obs,GSVA,how='left', left_index=True, right_index=True)
T2.obs=c
sc.pl.umap(T2, color=['Exhaustion','Cytotoxicity'], s=10, frameon=False, ncols=3, vmax='p99',save='featureplot.png')
cyto=T2[(T2.obs['celltype'] == 'GD')|(T2.obs['celltype'] == 'CD8_CTL'),:]
import re
r = re.compile('^CD8')
vmatch = np.vectorize(lambda x:bool(r.match(x)))
cyto=T2[(T2.obs['Cytotoxicity']>0)&(vmatch(T2.obs['celltype'])),:]
cyto=T2[(T2.obs['Exhaustion']>0)&(vmatch(T2.obs['celltype'])),:]
cyto=T2[(T2.obs['Exhaustion']>0),:]
cyto=T2[(T2.obs['Cytotoxicity']>0),:]
cyto=T2[(T2.obs['Cytotoxicity']>0),:]
cyto=T2[(T2.obs['Cytotoxicity']>0)&((T2.obs['celltype']=="GD")),:]
cyto=T2[(T2.obs['Cytotoxicity']>0)&((T2.obs['celltype']=="Unconv")),:]
cyto2=T2[(T2.obs['Cytotoxicity']>0)&((T2.obs['celltype']=="CD4_NOS")|(T2.obs['celltype']=="CD4_naive")|(T2.obs['celltype']=="Treg")),:]
pd.DataFrame([cyto.obs['Cytotoxicity'],cyto.obs['celltype']]).transpose().to_csv('/data/AMLscRNA/mergeblood/Tcell/check', index=False,header=True,sep="\t")

groupby_order = ["HD","RESb","RESa","UNb","UNa"]
dfa = sc.get.obs_df(cyto, ['Cytotoxicity','com'], use_raw=False)
df1 = dfa.melt(id_vars=["com"], value_vars="Cytotoxicity")
df1.rename(columns = {'com':'sample','value':'Cytotoxicity'}, inplace = True)
x = "sample"
y = "Cytotoxicity"
fig,ax = plt.subplots(figsize=(5,4),dpi=100,facecolor="w")
import seaborn as sns
ax = sns.boxplot(data=df1, x=x, y=y, order=groupby_order)
pairs=[("HD", "RESb"), ("HD", "UNb"), ("RESb","RESa"), ("UNb","UNa")]
from statannotations.Annotator import Annotator
annotator = Annotator(ax, pairs, data=df1, x=x, y=y, order=groupby_order)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.02,line_width=1)
annotator.apply_and_annotate()
ax.tick_params(which='major',direction='in',length=3,width=1.,labelsize=14,bottom=False)
for spine in ["top","left","right"]:
  ax.spines[spine].set_visible(False)

ax.spines['bottom'].set_linewidth(1)
ax.grid(axis='y',ls='--',c='gray')
ax.set_axisbelow(True)
import matplotlib.pyplot as plt
fig.tight_layout(pad=1, w_pad=1, h_pad=1)  
plt.savefig("Cytotoxicity_CD8.pdf")

T2.write_h5ad('T2_annotate.h5ad')

dfa = sc.get.obs_df(cyto, ['Exhaustion','com'], use_raw=False)
df1 = dfa.melt(id_vars=["com"], value_vars="Exhaustion")
df1.rename(columns = {'com':'CD8','value':'Exhaustion'}, inplace = True)
x = "CD8"
y = "Exhaustion"
fig,ax = plt.subplots(figsize=(5,4),dpi=100,facecolor="w")
import seaborn as sns
ax = sns.boxplot(data=df1, x=x, y=y, order=groupby_order)
pairs=[("HD", "RESb"), ("HD", "UNb"), ("RESb","RESa"), ("UNb","UNa")]
from statannotations.Annotator import Annotator
annotator = Annotator(ax, pairs, data=df1, x=x, y=y, order=groupby_order)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.02,line_width=1)
annotator.apply_and_annotate()
ax.tick_params(which='major',direction='in',length=3,width=1.,labelsize=14,bottom=False)
for spine in ["top","left","right"]:
  ax.spines[spine].set_visible(False)

ax.spines['bottom'].set_linewidth(1)
ax.grid(axis='y',ls='--',c='gray')
ax.set_axisbelow(True)
import matplotlib.pyplot as plt
fig.tight_layout(pad=1, w_pad=1, h_pad=1) 
plt.savefig("Exhaustion_CD8.pdf")
