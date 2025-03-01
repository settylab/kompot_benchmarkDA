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

def replicate_normalize_densities(sample_densities, replicate):
    sample_likelihoods = sample_densities.copy()
    for rep in list(replicate):
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        sample_likelihoods[curr_cols] = normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods


def runMELD(adata,k,sample_col, label_col, embedding_method,beta,meld_mode, dm_comp = 10):
    """
    meld_mode: if meld_mode == "optimal", using optimal values of meld. If the meld mode  == "default", dont specifically set the value for beta
    """
    # add sample and label dataframe to adata
    samplem = pd.DataFrame(index=pd.Series(adata.obs[sample_col]).unique())
    samplem.loc[:,label_col] = \
        adata.obs[[sample_col, label_col]].groupby(by=sample_col).aggregate(lambda x: x[0])
    adata.uns['samplem'] = samplem
    
    if embedding_method == "PCA":
        if meld_mode == "optimal":
            # run meld method
            if adata.n_vars <= 50:
                G = gt.Graph(adata.X, knn=k, use_pygsp=True)
            else:
                if 'X_pca' not in adata.obsm:
                    # perform pca and use it to generate graph
                    sc.tl.pca(adata, n_comps=50)
                G = gt.Graph(adata.obsm['X_pca_batch'], knn=k, use_pygsp=True)
            
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
        elif meld_mode == "default":
            # run meld method
            if adata.n_vars <= 50:
                G = gt.Graph(adata.X, knn=k, use_pygsp=True)
            else:
                if 'X_pca' not in adata.obsm:
                    # perform pca and use it to generate graph
                    sc.tl.pca(adata, n_comps=50)
                G = gt.Graph(adata.obsm['X_pca_batch'], knn=k, use_pygsp=True)
            
            meld_op = meld.MELD()
            # generate the densities of each sample
            sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[sample_col])
            # normalize the densities for each replicate
            replicates = samplem.index.map(lambda x: x.split('_')[-1]).unique()
            sample_likelihoods = replicate_normalize_densities(sample_densities, replicates)
            # average the likelihoods w.r.t conditions
            obj_cond = sorted(samplem[label_col].unique())[-1]
            obj_cond_columns = samplem.loc[samplem[label_col] == obj_cond,:].index.to_list()
            sample_likelihoods = sample_likelihoods[obj_cond_columns].mean(axis=1)
    
    elif embedding_method == "DiffusionMap":
            # run meld method
        if adata.n_vars <= 50:
            G = gt.Graph(adata.X, knn=k, use_pygsp=True)
        else:
            if "DM_EigenVectors" not in adata.obsm:
                palantir.utils.run_diffusion_maps(
                    adata, n_components=dm_comp, pca_key="X_pca_batch"
                )
            G = gt.Graph(adata.obsm["DM_EigenVectors"], knn=k, use_pygsp=True)
    
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

def calculate_outcome(adata,true_label,predicted_label):
    # calculate the metrics, TP, FP, TN, FN and so on.
    TP = ((adata.obs[true_label] == adata.obs[predicted_label]) & (adata.obs[predicted_label] != 'NotDA')).sum()
    FP = ((adata.obs[true_label] != adata.obs[predicted_label]) & (adata.obs[predicted_label] != 'NotDA')).sum()
    FN = ((adata.obs[true_label] != adata.obs[predicted_label]) & (adata.obs[predicted_label] == 'NotDA')).sum()
    TN = ((adata.obs[true_label] == adata.obs[predicted_label]) & (adata.obs[predicted_label] == 'NotDA')).sum()
    
    metrics_dic = {
        'TP': TP, 'FP': FP,
        'FN': FN, 'TN': TN,
        'TPR': [TP / (TP + FN)],
        'FPR': [FP / (FP + TN)],
        'TNR': [TN / (TN + FP)],
        'FNR': [FN / (FN + TP)],
        'FDR': [FP / (TP + FP)],
        'Precision': [TP / (TP + FP)],
        'Power': [1 - FN / (FN + TP)],
        'Accuracy': [(TP + TN) / (TP + TN + FP + FN)]
    }
    
    return pd.DataFrame(metrics_dic, index=['metric'])


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


def modified_auprc(adata, lfc):
    true_label = adata.obs["true_labels"]
    true_label_transformed_neg = (true_label == "NegLFC")
    true_label_transformed_pos = (true_label == "PosLFC")
    
    # Check if there are both positive and negative samples in the transformed arrays
    if np.any(true_label_transformed_neg) and np.any(~true_label_transformed_neg):
        auprc_score_neg = average_precision_score(true_label_transformed_neg, lfc)
    else:
        auprc_score_neg = np.nan
        
    if np.any(true_label_transformed_pos) and np.any(~true_label_transformed_pos):
        auprc_score_pos = average_precision_score(true_label_transformed_pos, lfc)
    else:
        auprc_score_pos = np.nan
    # except ValueError as e:
    #     print(f"ValueError: {e}")
    #     auprc_score_neg = np.nan
    #     auprc_score_pos = np.nan
    
    auprc_score = np.nanmean([auprc_score_neg, auprc_score_pos])
    return auprc_score,auprc_score_neg, auprc_score_pos

def modified_auroc(adata, lfc):
    true_label = adata.obs["true_labels"]
    true_label_transformed_neg = (true_label == "NegLFC")
    true_label_transformed_pos = (true_label == "PosLFC")
    

    if np.any(true_label_transformed_neg) and np.any(~true_label_transformed_neg):
        auc_score_neg = roc_auc_score(-true_label_transformed_neg, lfc)
        # if auc_score_neg <= 0.5:
        #     auc_score_neg = 1 - auc_score_neg
    else:
        auc_score_neg = np.nan

    if np.any(true_label_transformed_pos) and np.any(~true_label_transformed_pos):
        auc_score_pos = roc_auc_score(-true_label_transformed_pos, lfc)
        # if auc_score_pos <= 0.5:
        #     auc_score_pos = 1 - auc_score_pos
    else:
        auc_score_pos = np.nan
    # except ValueError as e:
    #     print(f"ValueError: {e}")
    #     auc_score_neg = np.nan
    #     auc_score_pos = np.nan
    
    auc_score = np.nanmean([auc_score_neg, auc_score_pos])
    return auc_score,auc_score_neg, auc_score_pos
    