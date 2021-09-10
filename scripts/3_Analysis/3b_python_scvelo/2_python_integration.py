## ANNDATA(PYTHON) INTEGRATION

import scvelo as scv
import pandas as pd
import anndata as ad

scv.logging.print_version()
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params(’scvelo’)

sdata = ad.read_mtx(’splicedMM.mtx’)
sdata.layers[’spliced’] = sdata.X
ndata = ad.read_mtx("unsplicedMM.mtx")
sdata.layers[’unspliced’] = ndata.X

pd_obs = pd.read_csv("seurat.obs.csv")
pd_obsm = pd.read_csv("seurat.obsm.csv")
sdata.obs = pd_obs
sdata.obsm["X_umap"] = pd_obsm.values

with open(’genenames.csv’, ’r’) as fil:
genes = fil.read().split(’\n’)
del genes[0]
del genes[-1]
genes
len(genes)
sdata.var_names = genes
