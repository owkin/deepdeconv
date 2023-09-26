"""scVI trainer debugger."""
# %%
import scanpy as sc
import scvi
import anndata as ad
import pandas as pd
import torch
import numpy

# %%
# load toy dataset
adata = scvi.data.heart_cell_atlas_subsampled()
sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="cell_source",
)

# %%
# Create and train model
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["cell_type"], #"cell_source", "donor"],
    # continuous_covariate_keys=["percent_mito", "percent_ribo"],
)
model = scvi.model.SCVI(adata)
model.view_anndata_setup()
model.train(max_epochs=10, batch_size=1024)

print("ok")
###
