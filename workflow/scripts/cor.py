
import numpy as np
import logging
import pandas as pd
import scipy.stats
import anndata

def weighted_pearson_correlation(A, B, wt, eps=1e-16):
        wt = np.ravel(wt) / np.sum(wt)
        Am = (A.T - np.dot(A, wt)).T
        Bm = (B.T - np.dot(B, wt)).T
        Am_sum_sq = np.dot(Am * Am, wt)
        Bm_sum_sq = np.dot(Bm * Bm, wt)
        numer = np.dot(Am * Bm, wt)
        denom = np.sqrt(Am_sum_sq * Bm_sum_sq)
        cor = np.divide(numer, denom, out=np.zeros_like(numer), where=denom != 0)
        return np.clip(cor, -1+eps, 1-eps)


def index_of(values, idx):
        sort_idx = np.argsort(idx)
        values_idx = sort_idx[np.searchsorted(idx, values, sorter=sort_idx)]
        assert np.all(idx[values_idx] == values)
        return values_idx

def correlate(rna, mark, interactions, weight, batch_size=10000, power=1, spearman=True, eps=1e-16):
        logger = logging.getLogger("epimap.linking")
        out = np.zeros(interactions.shape[0])
        g_idx = index_of(interactions["gene"].values, rna.var.index.values)
        e_idx = index_of(interactions["enh"].values, mark.var.index.values)
        for i_begin in range(0, interactions.shape[0], batch_size):
                i_end = min(i_begin + batch_size, interactions.shape[0])
                logger.info("Finding correlations %d -> %d" % (i_begin, i_end))
                cur_gene = g_idx[i_begin:i_end]
                cur_enh = e_idx[i_begin:i_end]
                R = rna.obsm["SVD_U"] @ np.diag(rna.uns["SVD_S"]**power) @ rna.varm["SVD_V"][cur_gene,:].T
                M = mark.obsm["SVD_U"] @ np.diag(mark.uns["SVD_S"]**power) @ mark.varm["SVD_V"][cur_enh,:].T
                if spearman:
                        R = scipy.stats.rankdata(R, axis=0)
                        M = scipy.stats.rankdata(M, axis=0)
                out[i_begin:i_end] = weighted_pearson_correlation(M.T, R.T, weight, eps=eps)
        return out

def run_correlation(rna, mark, pinteractions, ninteractions, weight=None, power=1. eps=1e-16, spearman=False, batch_size=10000):
        logger = logging.getLogger("epimap.linking")
        comm_samples = np.intersect1d(rna.obs.index, mark.obs.index)
        wt = np.ones(len(comm_samples))
        if weight is not None and weight in mark.obs.index:
                logger.info("Finding weight for BSSID %s" % weight)
                cov = mark[weight,:].obsm["SVD_U"] @ mark[comm_samples,:].obsm["SVD_U"].T
                Ncol = np.linalg.norm(mark[weight,:].obsm["SVD_U"], ord=2) + eps
                cov = np.ravel(np.asarray(cov)) / Ncol
                Nrow = np.linalg.norm(mark[comm_samples,:].obsm["SVD_U"], ord=2, axis=1) + eps
                wt = cov / Nrow
                wt = np.clip(np.asarray(wt), -1, 1)
                wt = wt ** 2
        M = mark[comm_samples,:]
        R = rna[comm_samples,:]
        logger.info("+ correlation")
        P = correlate(R, M, pinteractions, wt, batch_size=batch_size, spearman=spearman, power=power, eps=eps)
        logger.info("- correlation")
        N = correlate(R, M, ninteractions, wt, batch_size=batch_size, spearman=spearman, power=power, eps=eps)
        return (P, N)
