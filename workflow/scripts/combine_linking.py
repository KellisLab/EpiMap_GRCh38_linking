#!/usr/bin/env python3
import argparse
import pandas as pd


if __name__ == "__main__":
        ap = argparse.ArgumentParser()
        ap.add_argument("--bssid", required=True)
        ap.add_argument("--chromhmm", required=True)
        ap.add_argument("--linking", nargs="+")
        ap.add_argument("-o", "--output", required=True)
        args = vars(ap.parse_args())
        df = pd.concat([pd.read_csv(fname, header=None, sep="\t") for fname in args["linking"]])
        df.to_csv(args["output"], sep="\t", header=False, index=False)
