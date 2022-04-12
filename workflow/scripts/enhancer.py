import anndata
import numpy as np
import pandas as pd
import logging
import sys

def load_enh(chromhmm_h5ad, BSSID, state, min_overlap=0):
        logger = logging.getLogger("epimap.linking")
        chmm = anndata.read(chromhmm_h5ad)
        ## Subset ChromHMM matrix to active sample and state
        if BSSID not in chmm.layers:
                logger.critical("BSSID %s is not a valid biosample" % BSSID)
                sys.exit(1)
        ## now only filter for layer
        chmm = anndata.AnnData(chmm.layers[BSSID], obs=chmm.obs, var=chmm.var, dtype=chmm.layers[BSSID].dtype)
        if state not in chmm.obs.index.values:
                logger.critical("ChromHMM state %s must be one of: [%s]" % (state, ", ".join(chmm.obs.index.tolist())))
                sys.exit(1)
        pct = np.ravel(chmm[state,:].X.todense())
        enh = chmm.var.loc[pct >= min_overlap,:].copy()
        logger.info("Loaded %d enhancers" % (enh.shape[0]))
        enh["ChromHMM_overlap"] = np.ravel(chmm[state,enh.index.values].X.todense())
        enh["ChromHMM_state"] = state
        return enh
