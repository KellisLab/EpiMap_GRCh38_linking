#!/usr/bin/env python3
import anndata
import sys
import argparse 
import numpy as np
import os, sys
import pandas as pd

def run(input_list, output, sample_list, bed, col="mean"):
	bssid_list = pd.read_csv(sample_list, header=None)[0].values
	ccre = pd.read_csv(bed, sep="\t", header=None)
	ccre.columns = ["chrom", "begin", "end", "dhs", "enh", "comment"]
	ccre.index = ccre["enh"].values

	adata = anndata.AnnData(np.zeros((len(bssid_list),
					  ccre.shape[0]), dtype="float64"),
				obs=pd.DataFrame(index=bssid_list),
				var=ccre)
	adata.uns["bigwig_aggr"] = col
	bad = 0
	for bssid in bssid_list:
		print(bssid)
		F = [x for x in input_list if os.path.basename(x).startswith(bssid)]
		if len(F) != 1:
			print("For bssid", bssid, "input list is:", F)
		assert len(F) == 1
		try:
			df = pd.read_csv(F[0], sep="\t", header=None)
		except pd.errors.EmptyDataError:
			print("Bad file:", F[0], ", deleting")
			os.unlink(F[0])
			bad += 1
		if df.shape[1] == 6:
			df.columns = ["name", "size", "covered", "sum", "mean0", "mean"]
		elif df.shape[1] == 8:
			df.columns = ["name", "size", "covered", "sum", "mean0", "mean", "min", "max"]
		else:
			raise RuntimeError("%s does not have proper column format" % F[0])
		if col not in df.columns:
			raise RuntimeError("Column %s does not exist for file %s" % (col, F[0]))
		adata[bssid,df["name"].values].X = df[col].values
	if bad > 0:
		print("One or more samples were not present")
		sys.exit(1)
	print("writing to:", str(output))
	adata.write_h5ad(str(output))
	return adata

if __name__ == "__main__":
	if "snakemake" in globals():
		tbl = vars(snakemake)
		input_list = tbl["input"][2:]
		bed = tbl["input"][0]
		sample_list = tbl["input"][1]
		col = tbl["params"][0]
	else:
		ap = argparse.ArgumentParser()
		ap.add_argument("-i", "--input", nargs="+", required=True)
		ap.add_argument("-b", "--bed", required=True)
		ap.add_argument("-s", "--sample-list", required=True)
		ap.add_argument("-o", "--output", required=True)
		ap.add_argument("-c", "--column", default="mean")
		tbl = vars(ap.parse_args())
		input_list = tbl["input"]
		bed = tbl["bed"]
		sample_list = tbl["sample_list"]
		col = tbl["column"]
	run(input_list=input_list, output=tbl["output"], col=col,
	    sample_list=sample_list, bed=bed)
