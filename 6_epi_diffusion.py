#!/usr/bin/env python
# coding: utf-8

### DPT ###
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl

adata = sc.read_h5ad("epi.h5ad")
adata.obs['cell.type'].value_counts()

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15, random_state= 42)
sc.tl.diffmap(adata,n_comps=20)
sc.pl.diffmap(adata, color='cell.type', legend_loc='on data')
sc.pl.umap(adata, color=['day','cell.type'], legend_loc='on data')

adata.uns['iroot'] = np.flatnonzero(adata.obs['cell.type']  == 'Precursor')[0]
sc.tl.dpt(adata)
sc.pl.embedding(adata,'X_diffmap',color=['dpt_pseudotime','cell.type'])

pseudotime = adata.obs['dpt_pseudotime']
pseudotime.to_csv('pseudotime.csv',index=True)


