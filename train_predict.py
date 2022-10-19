import os
os.chdir('/data/AMLscRNA/mergeblood/20220501/undiff')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import scanpy as sc
import pandas as pd
import torch #pyro-ppl-1.6.0 torch-1.8.0
import scarches as sca 
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import scvi #0.8.1
condition_key = 'orig.ident'
cell_type_key = 'CellTypeN'
import os
os.chdir('/data/AMLscRNA/mergeblood/20220501/undiff')

vae_epochs = 500
scanvi_epochs = 200
surgery_epochs = 500
early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
early_stopping_kwargs_scanvi = {
    "early_stopping_metric": "accuracy",
    "save_best_state_metric": "accuracy",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
early_stopping_kwargs_surgery = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

ref = sc.read_loom("/data/yifan_AMLdata/dataset/Cell2019/raw/cell2019AMLfiltered.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
data = sc.read_h5ad("/data/AMLscRNA/mergeblood/20220501/undiff/undiff0501.h5ad")
adata_merge = ref.concatenate(data,batch_key='batch')  
adata_merge.obs['batch']=adata_merge.obs['batch'].cat.add_categories(['vanGalen','AML10']) 
adata_merge.obs.loc[adata_merge.obs['batch']=='0','batch']='vanGalen'
adata_merge.obs.loc[adata_merge.obs['batch']=='1','batch']='AML10'
adata_merge.obs['batch']=adata_merge.obs['batch'].cat.remove_unused_categories()
query=adata_merge[adata_merge.obs['CellTypeN'].isnull(),:]
ref=adata_merge[adata_merge.obs['CellTypeN'].notnull(),:]
sc.pp.log1p(ref)
sc.pp.highly_variable_genes(ref, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key="orig.ident")
var_genes_all = ref.var.highly_variable
print("Highly variable genes: %d"%sum(var_genes_all))
ref = ref[:, ref.var.highly_variable]
ref = sca.dataset.setup_anndata(ref, layer='counts', batch_key="orig.ident", labels_key='CellTypeN', copy=True)
vae = sca.models.SCANVI(
    ref,
    "Unknown",
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train(
    n_epochs_unsupervised=vae_epochs,
    n_epochs_semisupervised=scanvi_epochs,
    unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
    semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"],
                                       early_stopping_kwargs=early_stopping_kwargs_scanvi),
    frequency=1
)
reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["cell_type"] = ref.obs['CellTypeN'].tolist()
reference_latent.obs["batch"] = ref.obs['orig.ident'].tolist()
reference_latent.obs["Cohort"] = 'vanGalen'
sc.pp.neighbors(reference_latent, n_neighbors=10)
sc.tl.umap(reference_latent)
reference_latent.obs['predictions'] = vae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))
ref_path = '/data/AMLscRNA/mergeblood/20220501/undiff'
vae.save(ref_path, overwrite=True,save_anndata=True)
reference_latent.write_h5ad('reference0501.h5ad')

#predict
query= query.copy()
query.X = query.layers['counts']
HVGs = ref.var_names
query = query[:,HVGs]
query.obs[cell_type_key] = vae.unlabeled_category_

model = sca.models.SCANVI.load_query_data(query,ref_path,freeze_dropout = True)
model._unlabeled_indices = np.arange(query.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))
model.train(
    n_epochs_semisupervised=surgery_epochs,
    train_base_model=False,
    semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy", "elbo"],
                                       weight_decay=0,
                                       early_stopping_kwargs=early_stopping_kwargs_surgery
                                      ),
    frequency=1
)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['predictions'] = model.predict()
query_latent.obs['batch'] = query.obs['orig.ident'].tolist()

sc.pp.neighbors(query_latent)
sc.tl.umap(query_latent)
surgery_path = '/data/AMLscRNA/mergeblood/20220501/undiff/pre'
model.save(surgery_path, overwrite=True)

adata_full = ref.concatenate(query)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))

full_latent.obs['batch'] = adata_full.obs[condition_key].tolist()
full_latent.obs['state'] = adata_full.obs['state'].tolist()
full_latent.obs['Cohort'] = adata_full.obs['batch'].tolist()
full_latent.obs['predictions'] = model.predict(adata=adata_full)
full_latent.obs.index = adata_full.obs.index  

sc.pp.neighbors(full_latent, n_neighbors=30)
sc.tl.umap(full_latent, min_dist=0.2)
full_latent.write_h5ad('combined_ref_10p_0501.h5ad')
full_latent.obs[['Cohort','batch','state','predictions']].to_csv("scArches_PVG_DAV_CellType_0501.csv", index_label='Cell')


full_latent = sc.read_h5ad("/data/AMLscRNA/mergeblood/20220501/undiff/combined_ref_10p_0501.h5ad")
full_latent.obs['Cohort']=full_latent.obs['Cohort'].cat.add_categories(['vanGalen','AML10'])
full_latent.obs.loc[full_latent.obs['Cohort']=='0','Cohort']='vanGalen'
full_latent.obs.loc[full_latent.obs['Cohort']=='1','Cohort']='AML10'
full_latent.obs['Cohort']=full_latent.obs['Cohort'].cat.remove_unused_categories() 
from plotly.graph_objs import Scatter,Layout
import plotly
import plotly.offline as py
import plotly.graph_objs as go
import seaborn as sns
sns.set(rc={"figure.figsize":(14, 6)})
fig, (ax1, ax2) = plt.subplots(1, 2)
sc.pl.umap(full_latent[full_latent.obs.loc[full_latent.obs.Cohort.isin(['vanGalen'])].index.tolist(), ], color= ['predictions'], title=['vanGalen'], wspace = 0.6, ax=ax1,frameon=False,legend_loc="None",show=False,size=3)
sc.pl.umap(full_latent[full_latent.obs.loc[full_latent.obs.Cohort.isin(['AML10'])].index.tolist(), ], color= ['predictions'], title=['AML10'], wspace = 0.6, ax=ax2,frameon=False,legend_loc="right margin",show=False,size=3)
fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)  
plt.savefig("/data/AMLscRNA/mergeblood/20220501/undiff/umapsplit_vanGalen_DAV_new0502.pdf")
sc.pl.umap(full_latent, color= ['predictions'],  wspace = 0.6, ax=ax1,frameon=False,legend_loc="None",size=3)
