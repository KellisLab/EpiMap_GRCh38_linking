#!/usr/bin/env python3
import anndata
import sys
import argparse 
import numpy as np
import os
import pandas as pd
import pyranges
import requests
import gzip
from io import StringIO
import scipy.sparse

def download_chromhmm(bssid, url_base):
	r = requests.get(url_base + "/" + bssid + "_18_CALLS_segments.bed.gz")
	sio = StringIO(gzip.decompress(r.content).decode('utf-8'))
	df = pd.read_csv(sio,sep="\t", header=None)
	df.index = ["%s:%d-%d" % (x, y, z) for x, y, z in zip(df[0].values, df[1].values, df[2].values)]
	return pyranges.from_dict({"Chromosome": df[0].values,
				   "Start": df[1].values,
				   "End": df[2].values,
				   "State": df[3].values})

def index_of(values, idx):
	sort_idx = np.argsort(idx)
	values_idx = sort_idx[np.searchsorted(idx, values, sorter=sort_idx)]
	assert np.all(idx[values_idx] == values)
	return values_idx

	
def get_chmm_overlaps(bssid, url_base, ccre, adata):
	ccre_pr = pyranges.from_dict({"Chromosome": ccre["chrom"].values,
				      "Start": ccre["begin"].values,
				      "End": ccre["end"].values,
				      "CCRE": ccre["enh"].values})
	chmm_pr = download_chromhmm(bssid, url_base)
	join_df = ccre_pr.join(chmm_pr, report_overlap=True).df
	df = join_df.groupby(["CCRE", "State"]).agg(OverlapSum=("Overlap", "sum")).reset_index()
	df = pd.merge(df, df.groupby(["CCRE"]).agg(TotalSum=("OverlapSum", "sum")).reset_index())
	df["FracState"] = df["OverlapSum"] / df["TotalSum"]
	### finality: ChromHMM state X CCRE, layers=BSSID
	adata.layers[bssid] = np.zeros(adata.shape)
	ccre_idx = index_of(df["CCRE"].values, adata.var.index.values)
	state_idx = index_of(df["State"].values, adata.obs.index.values)
	adata.layers[bssid] = scipy.sparse.csr_matrix((df["FracState"].values, (state_idx, ccre_idx)),
						      shape=adata.shape)
	return adata


def run(input_list, output, chromhmm_url, states):
	ccre_df = None
	bssid_list = None
	for in_adata in input_list:
		adata = anndata.read(in_adata)
		if bssid_list is None:
			bssid_list = adata.obs.index.values
		else:
			assert np.all(bssid_list == adata.obs.index.values)
		if ccre_df is None:
			ccre_df = adata.var.copy()
		else:
			assert np.all(ccre_df.index.values == adata.var.index.values)
	adata = anndata.AnnData(obs=pd.DataFrame(index=states), var=ccre_df)
	for bssid in bssid_list:
		print(bssid)
		get_chmm_overlaps(bssid, chromhmm_url, ccre_df, adata)
	adata.write_h5ad(output)

if __name__ == "__main__":
	if "snakemake" in globals():
		tbl = vars(snakemake)
		chmm = tbl["params"][0]
		states = tbl["states"][1:]
		input_list = tbl["input"]
	else:
		ap = argparse.ArgumentParser()
		ap.add_argument("-i", "--input", nargs="+", required=True)
		ap.add_argument("-c", "--chromhmm", required=True)
		ap.add_argument("-o", "--output", required=True)
		ap.add_argument("-s", "--states", nargs="+", required=True)
		tbl = vars(ap.parse_args())
		input_list = tbl["input"]
		chmm = tbl["chromhmm"]
		states = tbl["states"]
	run(input_list=input_list, output=str(tbl["output"]), chromhmm_url=chmm, states=states)


