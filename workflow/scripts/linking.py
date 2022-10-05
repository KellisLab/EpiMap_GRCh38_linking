#!/usr/bin/env python3

from interaction import *
from cor import *
from expr import *
from mark import *
from enhancer import *

import scipy.sparse
import argparse
import numpy as np
import os, sys, re
import logging, json
from sklearn.linear_model import LogisticRegression
import ngboost
#from ngboost import NGBClassifier


def load(tss, rna, chromhmm, bssid, state, enh_overlap, enh_range, seed, mark_list, power=[], weight=True, square=True, neg_ratio=5., spearman=False, npc=833, transpose=False):
        logger = logging.getLogger("epimap.linking")
        ## now to enhancer
        enh = load_enh(chromhmm, BSSID=bssid, state=state, min_overlap=enh_overlap)

        rna = load_gene_info(tss, load_expr(rna, npc=npc, transpose=transpose))
        logger.info("Loaded %d genes with expression across %d samples" % (rna.shape[1], rna.shape[0]))

        pos_inter = positive_interactions(rna=rna.var, enh=enh, max_distance=enh_range)
        neg_inter = negative_interactions(enh=enh, pos_interactions=pos_inter, seed=seed, size=neg_ratio)

        # now to loading marks and correlation
        if power is None or len(power) == 0:
                power = [1]
        train_pos = np.zeros((pos_inter.shape[0], len(mark_list)*len(power)))
        train_neg = np.zeros((neg_inter.shape[0], len(mark_list)*len(power)))
        for i, mark_file in enumerate(mark_list):
                logger.info("Loading mark %s" % mark_file)
                mark = load_mark(mark_file, enh, npc=npc, transpose=transpose)
                if weight:
                        weight = bssid
                else:
                        weight = None
                for j, p in enumerate(power):
                        train_idx = i * len(power) + j
                        logger.info("Running correlation for power %.2f for mark %s" % (p, mark_file))
                        (P, N) = run_correlation(rna, mark, pos_inter, neg_inter, weight=weight, power=p, spearman=spearman)
                        train_pos[:, train_idx] = P
                        train_neg[:, train_idx] = N
                        del P, N
                del mark
        if square:
                train_pos = train_pos**2
                train_neg = train_neg**2
        train_pos = np.clip(train_pos, -1+1e-16, 1-1e-16)
        train_neg = np.clip(train_neg, -1+1e-16, 1-1e-16)
        train_pos = np.arctanh(train_pos)
        train_neg = np.arctanh(train_neg)
        return (pos_inter, neg_inter, train_pos, train_neg)

def train_logistic(dfp, dfn, pos, neg):
        """roadmap model"""
        dfp["score"] = 0
        pdbin = np.asarray(np.round(dfp["distance"].values / 1000), dtype="int")
        ndbin = np.asarray(np.round(dfn["distance"].values / 1000), dtype="int")
        updbin, pdbin_inv = np.unique(pdbin, return_inverse=True)
        undbin, ndbin_inv = np.unique(ndbin, return_inverse=True)
        for dbin in np.intersect1d(updbin, undbin):
                pi = np.isin(pdbin_inv, np.where(np.abs(dbin - updbin) < 5)[0])
                ni = np.isin(ndbin_inv, np.where(np.abs(dbin - undbin) < 5)[0])
                X = np.vstack((pos[pi,:], neg[ni,:]))
                y = np.hstack((np.ones(np.sum(pi)), np.zeros(np.sum(ni))))
                model = LogisticRegression(class_weight="balanced", warm_start=True).fit(X, y)
                pred_idx = pdbin_inv == np.where(dbin == updbin)[0]
                proba = model.predict_proba(pos[pred_idx,:])
                pos_class_idx = list(model.classes_).index(1)
                dfp.loc[pred_idx,"score"] = proba[:,pos_class_idx]
        return dfp

def abc_scale(dist, hic_power=0.87):
        offset = np.clip(np.abs(dist), 5000, np.Inf)
        scale = -4.80 + 11.63 * hic_power
        return np.exp(scale + -1 * hic_power * np.log(offset + 1))

def train_ngb(dfp, dfn, pos, neg, hic_power=0.87):
        y = np.hstack((np.ones(pos.shape[0]), np.zeros(neg.shape[0]))).astype(int)
        ### ABC weighting
        dfp["score"] = 0
        evalue = abc_scale(np.hstack((dfp["distance"].values, dfn["distance"].values)))
        model = ngboost.NGBClassifier(verbose=True)
        model.fit(np.hstack((np.vstack((pos, neg)),
                             evalue[:,None])),
                  y,
                  sample_weight=evalue)
        class_idx = 1
        evalue = abc_scale(dfp["distance"].values)
        dfp["score"] = model.predict_proba(np.hstack((pos, evalue[:,None])))[:,class_idx]
        return dfp

def dump(dfp, dfn, pos, neg, mark_list, output="output.h5ad"):
        out_base = os.path.join(os.path.dirname(output), os.path.basename(output).split('.')[0])
        out_pos = out_base + "_pos.h5ad"
        out_neg = out_base + "_neg.h5ad"
        dfp.index = ["P%d" % x for x in range(dfp.shape[0])]
        dfn.index = ["N%d" % x for x in range(dfn.shape[0])]
        mk = [os.path.basename(x).split(".")[0] for x in mark_list]
        anndata.AnnData(pos, obs=dfp, var=pd.DataFrame(index=mk)).write_h5ad(out_pos)
        anndata.AnnData(neg, obs=dfn, var=pd.DataFrame(index=mk)).write_h5ad(out_neg)
        return True

if __name__ == "__main__":
        ap = argparse.ArgumentParser(description="Generate gene-enhancer links")
        ap.add_argument("--rna", required=True, type=str, help="TSV matrix with id, gene, and log2fpkm columns")
        ap.add_argument("--mark", required=True, type=str, nargs="+", help="H5AD matrices containing histone mark/DNase data, (#BSSID x #DHS)")
        ap.add_argument("--bssid", required=True, type=str, help="Biosample ID, starts with 'BSS'")
        ap.add_argument("--state", required=True, type=str, help="ChromHMM state")
        ap.add_argument("--tss", required=True, type=str, help="BED6 file containing TSS info. Generated from gff3_to_bed.R")
        ap.add_argument("--chromhmm", required=True, type=str, help="H5AD file containing enhancer information including ChromHMM states and genomic locations")
        ap.add_argument("-o", "--output", required=True, dest="output", help="Output file (BED.GZ)")
        ap.add_argument("--cutoff", default=5/7, type=float, help="Cutoff of predicted positive class")
        ap.add_argument("--range", dest="enh_range", default=1000000, type=int, help="Distance from TSS at which enhancers can be called. Measured from mid-DHS")
        ap.add_argument("--overlap", dest="enh_overlap", default=0., type=float, help="Percent of enhancer-ChromHMM state overlap required for a potential enhancer")
        ap.add_argument("--neg-ratio", dest="neg_ratio", default=5., type=float, help="Percent of negative links to provide")
        ap.add_argument("--power", type=float, nargs="+")
        ap.add_argument("--weighted-cor", dest="weight", action="store_true")
        ap.add_argument("--no-weighted-cor", dest="weight", action="store_false")
        ap.add_argument("--dump", dest="dump", action="store_true")
        ap.add_argument("--model", dest="model", default="ngb", help="Machine learning model selection")
        ap.add_argument("--debug", dest="log_level", action="store_const", const=logging.DEBUG)
        ap.add_argument("--log-info", dest="log_level", action="store_const", const=logging.INFO)
        ap.add_argument("--quiet", dest="log_level", action="store_const", const=logging.WARNING)
        ap.add_argument("--seed", default=100, type=int, help="Random seed used in calculating negative E-G set")
        ap.add_argument("--npc", default=833, type=int, help="Number of PCs used in whitening and sample weighting")
        ap.add_argument("--spearman", dest="spearman", action="store_true", help="Use Spearman correlation instead of Pearson")
        ap.add_argument("--square", dest="square", action="store_true", help="Whether to use square correlations or just Pearson correlations")
        ap.add_argument("--no-square", dest="square", action="store_false", help="Whether to use square correlations or just Pearson correlations")
        ap.add_argument("--transpose", dest="transpose", action="store_true", help="Whether to use transposed PCA or just PCA")
        ap.add_argument("--no-transpose", dest="transpose", action="store_false", help="Whether to use transposed PCA or just PCA")
        ap.set_defaults(log_level=logging.INFO, dump=False, square=True, weight=False, spearman=False, transpose=False)
        args = vars(ap.parse_args())
        logging.basicConfig(level=args["log_level"], stream=sys.stdout)
        logger = logging.getLogger(args["output"])
        logger.info("Starting program")
        logger.info("Args: %s" % json.dumps(args))
        if args["dump"] and os.path.exists(args["output"]):
                logger.info("File %s exists already" % args["output"])
                sys.exit(1)
        dfp, dfn, pos, neg = load(tss=args["tss"], rna=args["rna"], chromhmm=args["chromhmm"], bssid=args["bssid"], state=args["state"], enh_overlap=args["enh_overlap"], enh_range=args["enh_range"], seed=args["seed"], mark_list=args["mark"], power=args["power"], weight=args["weight"], square=args["square"], spearman=args["spearman"], neg_ratio=args["neg_ratio"], npc=args["npc"], transpose=args["transpose"])
        if args["dump"]:
                dump(dfp, dfn, pos, neg, mark_list=args["mark"], output=args["output"])
        if args["model"] == "logistic":
                logger.info("Training logistic model")
                out = train_logistic(dfp, dfn, pos, neg)
        elif args["model"] in ["ngboost", "ngb"]:
                logger.info("Training NGBoost model")
                out = train_ngb(dfp, dfn, pos, neg)
        else:
                logger.critical("Model %s not implemented" % args["model"])
                sys.exit(1)
        cor_col = ["chrom", "begin", "end", "enh", "score", "distance", "gene", "gene_name", "tss", "ChromHMM_state"]
        out = out.loc[out["score"] >= args["cutoff"], cor_col]
        out.to_csv(args["output"],
                   index=False,
                   header=False,
                   sep="\t")
        logger.info("Done")
