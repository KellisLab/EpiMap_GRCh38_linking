import anndata
import numpy as np
import pandas as pd
from time import time
import logging
from sklearn.decomposition import PCA
import gc
import scipy.sparse

def load_mark(mark_h5, enh, transpose=False, npc=833):
        logger = logging.getLogger("epimap.linking")
        adata = anndata.read_h5ad(mark_h5, backed='r')
        adata = adata[:,enh.index.values] ## only load enhancers into memory that we want
        gc.collect()
        if "SVD_U" not in adata.obsm.keys():
                adata = adata.to_memory()
                if transpose:
                        U, s, VT = PCA()._fit(adata.X.T)
                        npc = min(len(s), npc)
                        adata.obsm["SVD_U"] = VT.T[:,range(npc)]
                        adata.uns["SVD_S"] = s[range(npc)]
                        adata.varm["SVD_V"] = U[:,range(npc)]
                else:
                        U, s, VT = PCA()._fit(adata.X)
                        npc = min(len(s), npc)
                        adata.obsm["SVD_U"] = U[:,range(npc)]
                        adata.uns["SVD_S"] = s[range(npc)]
                        adata.varm["SVD_V"] = VT.T[:,range(npc)]
        return anndata.AnnData(scipy.sparse.csr_matrix(([], ([], [])), shape=adata.shape),
                                obs=adata.obs,
                                var=adata.var,
                                obsm=adata.obsm,
                                uns=adata.uns,
                                varm=adata.varm)
