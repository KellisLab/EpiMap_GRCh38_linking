
import pyranges
import numpy as np
import pandas as pd
import logging
import sys

def positive_interactions(rna, enh, max_distance=1000000):
        logger = logging.getLogger("epimap.linking")
        r_gr = pyranges.from_dict({"Chromosome": rna["chrom"].values,
                                                           "Start": rna["tss"].values - max_distance,
                                                           "End": rna["tss"].values + max_distance,
                                                           "gene_index": rna.index.values })
        e_gr = pyranges.from_dict({"Chromosome": enh["chrom"].values,
                                                           "Start": enh["begin"].values,
                                                           "End": enh["end"].values,
                                                           "enh_index": enh.index.values })
        logger.info("Finding possible interactions between enhancers and genes")
        join_df = r_gr.join(e_gr).df
        if join_df.shape[0] == 0:
                logger.critical("Error: no possible E-G interactions with this ChromHMM state, biosample within %d base pairs" % max_distance)
                sys.exit(1)
        odf = enh.loc[join_df["enh_index"].values, :].copy()
        for col in ["gene", "tss", "strand", "gene_name"]:
                odf[col] = rna.loc[join_df["gene_index"].values,col].values
        odf["distance"] = (odf["begin"] + odf["end"])/2 - odf["tss"]
        odf["distance"] = odf["distance"] * (2 * (odf["strand"] == "+") - 1)
        odf = odf.loc[np.abs(odf["distance"]) <= max_distance,:]
        return odf

def negative_interactions(enh, pos_interactions, seed=100, col="enh", size=5):
        np.random.seed(seed)
        logger = logging.getLogger("epimap.linking")
        logger.info("Finding negative interactions")
        high = pos_interactions.shape[0]
        R = np.random.randint(low=0, high=high, size=int(round(high * size)))
        neg_interactions = pos_interactions.iloc[R,:]
        del R
        (uchrom, chrom_inv) = np.unique(neg_interactions["chrom"].values, return_inverse=True)
        for i, chrom in enumerate(uchrom):
                logger.info("Creating negative interactions on chrom %s" % chrom)
                echrom = enh.loc[enh["chrom"] != chrom,]
                I = np.random.randint(low=0, high=min(1, echrom.shape[0]), size=np.sum(i == chrom_inv))
                neg_interactions.loc[i == chrom_inv, col] = echrom[col].values[I]
        out_cols = np.union1d(["gene", "tss", "strand", "distance"], col)
        return neg_interactions[out_cols]
