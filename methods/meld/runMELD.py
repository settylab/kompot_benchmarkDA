import meld
import palantir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad

from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn import metrics

from sklearn.preprocessing import normalize
import graphtools as gt
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from lib import shared_embedding_utils
from lib.evaluation_metrics import calculate_outcome

def replicate_normalize_densities(sample_densities, replicate):
    sample_likelihoods = sample_densities.copy()
    for rep in list(replicate):
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        sample_likelihoods[curr_cols] = normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods


def runMELD(adata,k,sample_col, label_col, layer_embedding,beta, use_dm=False, dm_comp=10):
    # add sample and label dataframe to adata
    samplem = pd.DataFrame(index=pd.Series(adata.obs[sample_col]).unique())
    samplem.loc[:,label_col] = \
        adata.obs[[sample_col, label_col]].groupby(by=sample_col).aggregate(lambda x: x.iloc[0])
    adata.uns['samplem'] = samplem

    if adata.n_vars <= 50:
        G = gt.Graph(adata.X, knn=k, use_pygsp=True)
    else:
        # Use shared embedding utilities for consistency
        adata = shared_embedding_utils.ensure_batch_corrected_embeddings(
            adata, layer_embedding=layer_embedding, dm_comp=dm_comp
        )

        # Get the appropriate embedding matrix
        X, embedding_name = shared_embedding_utils.get_embedding_for_method(
            adata, use_dm=use_dm, dm_comp=dm_comp
        )

        G = gt.Graph(X, knn=k, use_pygsp=True)
        
    meld_op = meld.MELD(beta=beta)
        # generate the densities of each sample
    sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[sample_col])
        # normalize the densities for each replicate
    replicates = samplem.index.map(lambda x: x.split('_')[-1]).unique()
    sample_likelihoods = replicate_normalize_densities(sample_densities, replicates)
        # average the likelihoods w.r.t conditions
    obj_cond = sorted(samplem[label_col].unique())[-1]
    obj_cond_columns = samplem.loc[samplem[label_col] == obj_cond,:].index.to_list()
    sample_likelihoods = sample_likelihoods[obj_cond_columns].mean(axis=1)

    return sample_likelihoods.values,samplem


def threshold_meld(adata,sample_likelihoods):

    lower = sample_likelihoods.min() + 1e-8
    upper = sample_likelihoods.max() - 1e-8
    assert (lower < upper), \
        "lower bound is not less than upper bound"
    threshold_meld_res = np.arange(lower, upper, (upper - lower) / 100)
    
    return threshold_meld_res


def meld2output(meld_res, out_type="continuous", thresholds=None):
    if out_type == "continuous":
        da_cell = meld_res
    else:
        def get_da_cell(thres):
            isPos = meld_res > thres
            isNeg = meld_res < 1 - thres
            
            da = np.array(["NotDA"] * len(meld_res), dtype=object)
            da[isPos] = "PosLFC"
            da[isNeg] = "NegLFC"
            return da
        
        if isinstance(thresholds, float):
            da_cell = get_da_cell(thresholds)
        elif isinstance(thresholds, list) or isinstance(thresholds, np.ndarray):
            da_cell = []
            for _, thres in enumerate(thresholds):
                da_cell.append(get_da_cell(thres))
        else:
            raise RuntimeError("param: alphas can only support list or float")
    
    return da_cell


def get_performance_df_meld(adata,true_label_col,predicted_labels_meld,threshold_meld_res):
    new_adata = adata.copy()
    predicted_labels = {}
    performance_df = []
    for i in range(len(predicted_labels_meld)):  
        predicted_labels["predicted_label_meld_" + str(i)] = predicted_labels_meld[i]
    
    predicted_labels_df = pd.DataFrame(predicted_labels,index = new_adata.obs_names)

    new_adata.obs = pd.concat([new_adata.obs, predicted_labels_df],axis = 1)
    
    for j in range(len(predicted_labels_meld)):
        temp_performance = calculate_outcome(new_adata,true_label_col,"predicted_label_meld_"+str(j))
        performance_df.append(temp_performance)
    result_df = pd.concat(performance_df)
    result_df["threshold"] = list(threshold_meld_res)
    return result_df



def evaluation(result_performance):
    auc = metrics.auc(result_performance['FPR'].values, result_performance['TPR'].values)
    prc = metrics.auc(result_performance["TPR"].values,result_performance["Precision"].values)
    print("AUC: ", auc)
    print("PRC: ", prc)
    return auc,prc
