
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import logging
import re
import anndata
from sklearn.decomposition import PCA
### Wrapper around RNA-seq

def load_expr(tsv, transpose=False, npc=833):
        logger = logging.getLogger("epimap.linking")
        logger.info("Loading expression from %s" % tsv)
        tbl = pd.read_csv(tsv, sep="\t")
        (BSSID, I) = np.unique(tbl['id'].values, return_inverse=True)
        (GENES, G) = np.unique(tbl['gene'].values, return_inverse=True)
        mat = csr_matrix((tbl['log2fpkm'].values, (I, G)), shape=(BSSID.shape[0], GENES.shape[0])).todense()
        npc = min(npc, len(BSSID))
        logger.debug("Loaded expression matrix of size (%d,%d)" % mat.shape)
        adata = anndata.AnnData(np.asarray(mat), obs=pd.DataFrame(index=BSSID), var=pd.DataFrame(index=GENES), dtype=mat.dtype)
        if transpose:
                U, s, VT = PCA()._fit(mat.T)
                adata.obsm["SVD_U"] = VT.T[:,range(npc)]
                adata.uns["SVD_S"] = s[range(npc)]
                adata.varm["SVD_V"] = U[:,range(npc)]
        else:
                U, s, VT = PCA()._fit(mat)
                adata.obsm["SVD_U"] = U[:,range(npc)]
                adata.uns["SVD_S"] = s[range(npc)]
                adata.varm["SVD_V"] = VT.T[:,range(npc)]
        adata.X = csr_matrix(([], ([], [])), shape=adata.shape) ## use PCA
        return adata

def load_gene_info(gene_info_bed, RNA_matrix):
        logger = logging.getLogger("epimap.linking")
        logger.info("Loading gene info from %s" % gene_info_bed)
        ginfo = pd.read_csv(gene_info_bed, sep="\t", header=None)
        logger.debug("Loaded %d gene TSS" % ginfo.shape[0])
        ginfo = ginfo.iloc[:,0:6]
        ginfo.columns = ["chrom", "begin", "end", "gene_name", "gene_id", "strand"]
        ginfo["gene_name"] = [re.sub("[.][_0-9]+$", "", str(g)) for g in ginfo["gene_name"].values]
        ginfo["gene_id"]   = [re.sub("[.][_0-9]+$", "", g) for g in ginfo["gene_id"].values]
        L_gn = len(np.intersect1d(ginfo["gene_name"].values, RNA_matrix.var.index.values))
        L_gi = len(np.intersect1d(ginfo["gene_id"].values, RNA_matrix.var.index.values))
        ### Decide if matrix is ENSEMBL or SYMBOL
        if L_gn > L_gi:
                ginfo["gene"] = ginfo["gene_name"].values
        else:
                ginfo["gene"] = ginfo["gene_id"].values
        ginfo = ginfo.drop_duplicates("gene")
        ginfo.index = ginfo["gene"].values
        logger.debug("Found TSS for %d unique genes" % ginfo.shape[0])
        ginfo["tss"] = ginfo["begin"].values
        ginfo.loc[ginfo["strand"] != "+", "tss"] = ginfo["end"]
        assert np.sum(np.diff(ginfo.loc[ginfo["strand"] == "+", ["tss", "begin"]].values, axis=1)) == 0
        assert np.sum(np.diff(ginfo.loc[ginfo["strand"] != "+", ["tss", "end"]].values, axis=1)) == 0
        ### inner join on indices
        glist = np.intersect1d(ginfo.index.values, RNA_matrix.var.index.values)
        RNA_matrix = RNA_matrix[:,glist].copy()
        RNA_matrix.var = ginfo.loc[glist,:]
        return RNA_matrix
