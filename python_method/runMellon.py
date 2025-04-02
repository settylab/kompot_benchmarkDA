import mellon
import palantir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad

from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn import metrics

import logging
logger = logging.getLogger("mellon")

def compute_z_score(dens1, dens2, var1, var2, eps=1e-16):
    log_fold_change_mean = dens2 - dens1
    zscores = log_fold_change_mean / np.sqrt(var1 + var2 + eps)
    return log_fold_change_mean, zscores


def runMELLON(
    adata, mellon_d_method: str, norm_density:str,  label_col: str, ls_mode:str,layer_embedding,dm_comp: int = 10,  ls_factor: float = 1):

    conditions = adata.obs[label_col].unique()

    if adata.n_vars <= 50:
        if ls_mode == "PCA":
            if not isinstance(adata.X, np.ndarray):
                adata.obsm[f"{layer_embedding}_batch"] = adata.X.toarray()
            else:
                adata.obsm[f"{layer_embedding}_batch"] = adata.X
        elif ls_mode == "DM":
            adata.obsm[f"{layer_embedding}_batch"] = adata.obsm[f"{layer_embedding}_batch"]
    else:
        if f"{layer_embedding}_batch" not in adata.obsm:
            logging.info("Running PCA")
            # perform pca and use it to generate graph
            sc.tl.pca(adata, n_comps=50)
            adata.obsm[f"{layer_embedding}_batch"] = adata.obsm["X_pca"]
    if dm_comp > 0:
       if ("DM_EigenVectors" not in adata.obsm
            or adata.obsm["DM_EigenVectors"].shape[1] != dm_comp):
            palantir.utils.run_diffusion_maps(adata, n_components=dm_comp, pca_key=f"{layer_embedding}_batch")
            X = adata.obsm["DM_EigenVectors"]
            cov_func_curry = mellon.cov.Matern52
    else:
        X = adata.obsm[f"{layer_embedding}_batch"]
        if not isinstance(X, np.ndarray):
            X = X.to_numpy() 
        if ls_mode == "DM":
            ls_factor = ls_factor
            logger.info(f"ls_factor large, not change,ls_factor={ls_factor}")
        elif ls_mode == "PCA":
            ls_factor *= 2
            logger.info(f"ls_factor small, ls_factor={ls_factor}")
        cov_func_curry = mellon.cov.Matern52
    # X = adata.obsm["X_pca"]
    # X = adata.obsm["X_pca"]

    densities = list()
    variances = list()
    densities_norm = []
    d_list = []
    ls_list = []
    mu_list = []
    for condition in sorted(conditions):
        idx = adata.obs[label_col] == condition
        n = np.sum(idx)
        sub_X = X[idx, :]

        estimator = mellon.DensityEstimator(
            ls_factor=ls_factor,
            d_method = mellon_d_method,
            cov_func_curry=cov_func_curry,
            optimizer="advi",
            predictor_with_uncertainty=True
        )
        predictor = estimator.fit(sub_X).predict
        dens = predictor(X)
        dens_norm = predictor(X, normalize=True)
        # dens -= d*np.log(n) # normalization
        densities.append(dens)
        var = predictor.uncertainty(X)
        variances.append(var)
        densities_norm.append(dens_norm)
        d = estimator.d
        ls = estimator.ls
        mu = estimator.mu
    
        d_list.append(d)
        ls_list.append(ls)
        mu_list.append(mu)
        
    
    if norm_density == "Yes" and mellon_d_method == "fractal":
        log_fold_change_mean, zscores = compute_z_score(*densities_norm, *variances)
    elif norm_density == "No" and mellon_d_method == "fractal":
        log_fold_change_mean, zscores = compute_z_score(*densities, *variances)
    else:
        log_fold_change_mean, zscores = compute_z_score(*densities, *variances)

    return log_fold_change_mean, zscores



def runMELLON_synchronized(
    adata, mellon_d_method: str, norm_density:str, corrected:str, label_col: str, dm_comp,  ls_factor: float , ls_mode:str,layer_embedding):
    
    """
    This function is different from the above runMellon function, since it estimated hyper-parameters including mu for Gaussian Distribution, landscale factor and d.
    ls_mode: if PCA: ls small; if DM: ls large
    (Note: the mode of PCA or DM needs to be set manually, because for avoiding additional error, the DM eigenvectors replaced the PCA layer in anndata)
    """

    conditions = adata.obs[label_col].unique()

    if adata.n_vars <= 50:
        if ls_mode == "PCA":
            if not isinstance(adata.X, np.ndarray):
                adata.obsm[f"{layer_embedding}_batch"] = adata.X.toarray()
            else:
                adata.obsm[f"{layer_embedding}_batch"] = adata.X
        elif ls_mode == "DM":
            adata.obsm[f"{layer_embedding}_batch"] = adata.obsm[f"{layer_embedding}_batch"]
    else:
        if f"{layer_embedding}_batch" not in adata.obsm:
            logging.info("Running PCA")
            # perform pca and use it to generate graph
            sc.tl.pca(adata, n_comps=50)
            adata.obsm[f"{layer_embedding}_batch"] = adata.obsm["X_pca"]
            
    if dm_comp > 0:
       if ("DM_EigenVectors" not in adata.obsm
            or adata.obsm["DM_EigenVectors"].shape[1] != dm_comp):
            palantir.utils.run_diffusion_maps(adata, n_components=dm_comp, pca_key=f"{layer_embedding}_batch")
            X = adata.obsm["DM_EigenVectors"]
            cov_func_curry = mellon.cov.Matern52

    else:
        X = adata.obsm[f"{layer_embedding}_batch"]
        if not isinstance(X, np.ndarray):
            X = X.to_numpy() 
        if ls_mode == "DM":
            ls_factor = ls_factor
            logger.info(f"ls_factor large, not change,ls_factor={ls_factor}")
        elif ls_mode == "PCA":
            ls_factor *= 2
            logger.info(f"ls_factor small, ls_factor={ls_factor}")
        cov_func_curry = mellon.cov.Matern52
    # X = adata.obsm["X_pca"]
    
    # # compute hyper parameters
    # Compute hyperparameters
    d = mellon.parameters.compute_d_factal(X)
    logger.info(f"Computed d={d}")
    nn_distances = mellon.parameters.compute_nn_distances(X)
    logger.info("Computed nearest neighbor distances")
    base_ls = mellon.parameters.compute_ls(nn_distances)
    logger.info(f"Computed base_ls={base_ls}")
    if ls_mode == "DM":
        ls = ls_factor * base_ls
        logger.info(f"Computed DM layer ls={ls}")
    elif ls_mode == "PCA":
        ls = ls_factor * base_ls * 2 ** (1 / d)
        logger.info(f"Computed PCA layer ls={ls}")
    mu = mellon.parameters.compute_mu(nn_distances, d) - 5
    logger.info(f"Computed mu={mu}")
    landmarks_compute = mellon.parameters.compute_landmarks(X, n_landmarks=5_000)
    logger.info("Computed landmarks")

    #lk_num = mellon.parameters.compute_n_landmarks("sparse_cholesky", adata.n_obs,None)
    #lk = mellon.parameters.compute_landmarks(X,n_landmarks = 0)
    
    estimators = []
    for condition in sorted(conditions):
        logger.info(f"Processing condition {condition}.")
        idx = adata.obs[label_col] == condition
        sub_X = X[idx, :]
        estimator = mellon.DensityEstimator(
            ls=ls,
            cov_func_curry=cov_func_curry,
            landmarks=landmarks_compute,
            optimizer="advi",
            predictor_with_uncertainty=True,
            d=d,
            mu=mu,
        )
        estimator.fit(sub_X)
        estimators.append(estimator)

    logger.info("Computing densities")
    densities = [e.predict(X) for e in estimators]
    logger.info("Computing normalized densities")
    normalized_densities = [e.predict(X, normalize=True) for e in estimators]
    logger.info("Computing uncertainty")
    variances = [e.predict.uncertainty(X) for e in estimators]
        
    # var = predictor.uncertainty(X)
    # variances.append(var)

    # condition2_dens = densities[1]
    # condition1_dens = densities[0]
    # condition2_dens_norm = densities_norm[1]
    # condition1_dens_norm = densities_norm[0]




    if norm_density == "Yes" and mellon_d_method == "fractal":
        logger.info("Computing normalized log-fold change.")
        log_fold_change_mean, zscores = compute_z_score(*normalized_densities[:2], *variances[:2])
    elif norm_density == "No" and mellon_d_method == "fractal":
        logger.info("Computing log-fold change.")
        log_fold_change_mean, zscores = compute_z_score(*densities[:2], *variances[:2])
    else:
        logger.info("Computing log-fold change with default.")
        log_fold_change_mean, zscores = compute_z_score(*densities[:2], *variances[:2])

    
    if corrected == "Yes":
        logger.info("Applying corrective term")
        ddens = sum(densities) * d
        correction = np.exp(ddens - np.mean(ddens))
        log_fold_change_mean = log_fold_change_mean * correction
    else:
        if corrected == "No":
            logger.info("Applying uncorrective term, the log fold change mean not change")
            log_fold_change_mean = log_fold_change_mean
    
    return log_fold_change_mean, zscores
    #return condition1_dens,condition2_dens,log_fold_change_mean, zscores


def threshold_mellon(zscores):
    # get da cells with different thresholds
    lower = 0
    upper =np.max(np.abs(zscores)) - 1e-8
    #np.max(np.abs(zscores)) - 1e-8
    assert lower < upper, "lower bound is not less than upper bound"
    thresholds = np.arange(lower, upper, (upper - lower) / 100)
    return thresholds

def mellon2output(lfcs, zscores, out_type="continuous", thresholds=None):
    if out_type == "continuous":
        da_cell = zscores
    else:

        def get_da_cell(thres):
            issign = np.abs(zscores) > thres
            isPos = np.logical_and(issign, lfcs > 0)
            isNeg = np.logical_and(issign, lfcs < 0)

            da = np.array(["NotDA"] * len(zscores), dtype=object)
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


def get_performance_df(adata,true_label_col,predicted_labels_mellon,thresholds):
    new_adata = adata.copy()
    predicted_labels = {}
    performance_df = []
    for i in range(len(predicted_labels_mellon)):  
        predicted_labels["predicted_label_" + str(i)] = predicted_labels_mellon[i]
    
    predicted_labels_df = pd.DataFrame(predicted_labels,index = new_adata.obs_names)

    new_adata.obs = pd.concat([new_adata.obs, predicted_labels_df],axis = 1)
    
    for j in range(len(predicted_labels_mellon)):
        temp_performance = calculate_outcome(new_adata,true_label_col,"predicted_label_"+str(j))
        performance_df.append(temp_performance)
    result_df = pd.concat(performance_df)
    result_df["threshold"] = list(thresholds)
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
    
