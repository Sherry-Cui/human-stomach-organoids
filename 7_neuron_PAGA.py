# Extended Data Figure6
# import module
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.pyplot import rc_context
from matplotlib.backends.backend_pdf import PdfPages
 

# load file
fn = "Neuron.h5ad"
te = sc.read_h5ad(fn)

custom_palette = {
        'Migratory ENCC':"#6BD76BFF",
        'NE':'#A6CEE3',
        'NPC':"#1F78B4",
        'Neuron1':"#e7c23e",
        'Neuron2':"#FF4500"  ,
        'Premigratory ENCC': "#802268FF"}
sc.pl.umap(te, color=['cell.type'],palette=custom_palette)
sc.pp.neighbors(te,use_rep='X_umap',n_neighbors=100)
sc.tl.paga(te, groups='cell.type')
sc.pl.paga(te,colors='cell.type')

np.savetxt(fname="PAGA.csv",
           X=te.uns['paga']['connectivities'].todense(), delimiter=",")

te.write('PAGA.h5ad', compression='gzip')

te.obs['cell.type'].cat.categories
sc.pl.paga(te, threshold=0.086, show=False)
te.uns['paga']['connectivities'].todense()


sc.tl.draw_graph(te)
sc.pl.draw_graph(te,color=['cell.type'], legend_loc='on data')
with rc_context({'figure.figsize': (10, 25)}):
    sc.pl.paga_compare(te,basis='X_umap' ,threshold=0.4,title='', right_margin=0, size=10, 
                   edge_width_scale=1,min_edge_width=5,legend_fontsize=20, fontsize=10, frameon=False, edges=True,node_size_scale=20,plot=True,save='PAGA.1.pdf')


fig, axs = plt.subplots(figsize=(7, 5.5))
# point for each fig show=False
sc.pl.embedding(te,basis='X_umap',color='cell.type',ax=axs,show=False,size=30)
sc.pl.paga(te,color='cell.type',layout='fa',threshold=0.4,
           node_size_power=1,max_edge_width=2,min_edge_width=0.5,
           add_pos=False,pos=te.uns['paga']['pos'],fontsize=8,node_size_scale=0.05,
          edge_width_scale=0.5,ax=axs,plot=True,save='PAGA.pdf')

# SAVE PLOT
with plt.rc_context({'figure.figsize':(6.8,5),'figure.dpi':300}):
    with rc_context({'figure.figsize': (10, 25)}):
       sc.pl.paga_compare(te,basis='X_umap' ,threshold=0.086,title='', right_margin=0, size=10, 
                       edge_width_scale=1,min_edge_width=5,legend_fontsize=20, fontsize=10, frameon=False, edges=True,node_size_scale=20)

    plt.subplots_adjust(right=0.8)
    plt.savefig('PAGA.pdf',bbox_inches = 'tight')
