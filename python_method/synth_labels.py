
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r, ListVector

#import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

import pandas as pd
from rpy2.robjects import pandas2ri


import os
import sys



def cap_probabilities(adata,cond_probability,conditions, balance = "Yes",cap_enr=None):
    
    """
    Genrate condition 1 probability and condition 2 probability
    """
    
    if cap_enr is not None:
        cond_probability = np.where(cond_probability > cap_enr, cap_enr, cond_probability)
    
    if balance == "Yes":
        norm_factor = 0.5 * adata.n_obs / np.sum(cond_probability)
    elif balance == "No":
        norm_factor = 1
    condition1_prob = norm_factor * cond_probability
    condition2_prob = 1 - condition1_prob
    
    cond_probability_df = pd.DataFrame({
    conditions[0]: condition1_prob,
    conditions[1]: condition2_prob
    })
    #conditions,cond_probability_df,w_logit,centroid_distance
    return conditions,cond_probability_df


def label_condition_and_rep_labels(adata,cond_probability,seed):
           
    """
    Generated synth_labels based in the condition probability
    Based on n_batches to label batches
    Based on n_replicates and synth_labels to label samples.
    
    """
    # cond_probability = pd.DataFrame()
    # cond_probability.loc[:,0] = adata.obs["condition1_prob"]
    # cond_probability.loc[:,1] = adata.obs["condition2_prob"]
    # cond_probability.columns = conditions
    #np.random.seed(seed)
    # Suppress pandas2ri deprecation warning and convert DataFrame
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        # Convert the pandas DataFrame to an R data frame within a conversion context
        with localconverter(robjects.default_converter + pandas2ri.converter):
            r_cond_probability = robjects.conversion.py2rpy(cond_probability)

    # Assign the R data frame to an R variable in the global environment
    robjects.globalenv['cond_probability'] = r_cond_probability

    robjects.globalenv['seed'] = seed 
    # Execute the R code
    robjects.r('''
    set.seed(seed)
    synth_labels <- sapply(1:nrow(cond_probability), function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
    ''')

    # Convert the R vector back to a pandas Series
    with localconverter(robjects.default_converter + pandas2ri.converter):
        synth_labels = list(robjects.globalenv['synth_labels'])
    
#     replicates = [f"R{i}" for i in range(1, n_replicates + 1)]
    
#     batch_labels = [f"B{i}" for i in range(1, n_batches + 1) for _ in range(n_replicates)]
#     np.random.shuffle(batch_labels)
#     batches = batch_labels
    
#     synth_samples = [f"{label}_{rep}" for label, rep in zip(synth_labels, replicates * len(synth_labels))]
    
#     if n_batches > 1:
#     # Use a dictionary comprehension to map unique synth_samples to sorted unique labels
#         batch_labels_dict = {sample:batches[i] for i, sample in enumerate(sorted(set(synth_samples)))}
#     else:
#         # If there is only one batch, map all unique synth_samples to "B1"
#         batch_labels_dict = {sample: "B1" for sample in set(synth_samples)}
#     synth_batches = [batch_labels_dict[sample] for sample in synth_samples]
    
#     synth_batches_df = pd.DataFrame(synth_batches, index = synth_samples)
    
    adata.obs["synth_labels"] = synth_labels
    # adata.obs["synth_samples"] = synth_samples
    # if synth_samples == list(synth_batches_df.index):
    #     print("yes")
    # adata.obs["synth_batches"] = list(synth_batches_df.iloc[:,0])
    adata.obs["Condition1_prob"] = list(cond_probability.iloc[:,0])
    adata.obs["Condition2_prob"] = list(cond_probability.iloc[:,1])
    
    return adata


def label_condition_and_rep_other(adata,n_replicates, n_batches,seed):
           
    """
    Generated synth_labels based in the condition probability
    Based on n_batches to label batches
    Based on n_replicates and synth_labels to label samples.
    
    """
    # cond_probability = pd.DataFrame()
    # cond_probability.loc[:,0] = adata.obs["condition1_prob"]
    # cond_probability.loc[:,1] = adata.obs["condition2_prob"]
    np.random.seed(seed)
    
    # synth_labels = [
    #         np.random.choice(a=conditions,replace=False,p=cond_probability.iloc[i, :])
    #         for i in range(len(cond_probability))
    # ]
    
    synth_labels = adata.obs["synth_labels"]
    
    replicates = [f"R{i}" for i in range(1, n_replicates + 1)]
    
    batch_labels = [f"B{i}" for i in range(1, n_batches + 1) for _ in range(n_replicates)]
    np.random.shuffle(batch_labels)
    batches = batch_labels
    
    synth_samples = [f"{label}_{rep}" for label, rep in zip(synth_labels, replicates * len(synth_labels))]
    
    if n_batches > 1:
    # Use a dictionary comprehension to map unique synth_samples to sorted unique labels
        batch_labels_dict = {sample:batches[i] for i, sample in enumerate(sorted(set(synth_samples)))}
    else:
        # If there is only one batch, map all unique synth_samples to "B1"
        batch_labels_dict = {sample: "B1" for sample in set(synth_samples)}
    synth_batches = [batch_labels_dict[sample] for sample in synth_samples]
    
    synth_batches_df = pd.DataFrame(synth_batches, index = synth_samples)
    
    #adata.obs["synth_labels"] = synth_labels
    adata.obs["synth_samples"] = synth_samples
    if synth_samples == list(synth_batches_df.index):
        print("yes")
    adata.obs["synth_batches"] = list(synth_batches_df.iloc[:,0])

    
    return adata




def quantile_assign_label(adata,pop_col,pop_enr,pop):
    
    """
    Assign "true_labels" to adata, allow to perform mellon or meld and do the performance evaluation.
    """
    pop_tbl = adata.obs[pop_col].value_counts()
    
    total_cells = pop_tbl.sum()
    n_enr_cells = pop_tbl[pop]
    neutral_prop = (total_cells - n_enr_cells) / (total_cells*2)
    print(neutral_prop)
    # Determine the quantile based on the condition
    if pop_enr < 0.5:
        da_lower = np.quantile(adata.obs['Condition1_prob'], pop_tbl[pop_tbl.index == pop] / total_cells)[0]
        da_upper = (2*neutral_prop) - da_lower
    else:
        da_upper = np.quantile(adata.obs['Condition1_prob'], 1 - pop_tbl[pop_tbl.index == pop] / total_cells)[0]
        print(da_upper)

        da_lower = (2*neutral_prop) - da_upper
        print(da_lower)
        
    assert da_upper > da_lower, "da_upper must be greater than da_lower"
    
    bins = [-float('inf'), da_lower, da_upper, float('inf')]
    labels = ["PosLFC", "NotDA", "NegLFC"]
    adata.obs['true_labels'] = pd.cut(adata.obs['Condition1_prob'], bins, labels=labels)
    
    return adata
    

def quantile_assign_label_old(adata,pop_col,pop_enr,pop):
    
    """
    Assign "true_labels" to adata, allow to perform mellon or meld and do the performance evaluation.
    """
    pop_tbl = pd.DataFrame(adata.obs[pop_col].value_counts())
    
    # Determine the quantile based on the condition
    if pop_enr < 0.5:
        da_lower = np.quantile(adata.obs['Condition1_prob'], pop_tbl[pop_tbl.index == pop] / sum(pop_tbl.iloc[:,0]))
        da_upper = 1 - da_lower
    else:
        da_upper = np.quantile(adata.obs['Condition1_prob'], 1 - pop_tbl[pop_tbl.index == pop] / sum(pop_tbl.iloc[:,0]))
        da_lower = 1 - da_upper
        
    assert da_upper > da_lower, "da_upper must be greater than da_lower"
    
        # Assign the synthetic labels based on 'Condition2_prob'
    true_label = np.where(np.array(adata.obs['Condition2_prob']) < da_lower, 'NegLFC',
                                     np.where(np.array(adata.obs['Condition2_prob']) > da_upper, 'PosLFC', 'NotDA'))

    
    adata.obs['true_labels'] = true_label[0]
    # Replace spaces with underscores in pop if necessary
    pop = pop.replace(" ", "_")
    
    return adata
    
def add_batch_effect_pca(adata, layer_embedding, batch_col="synth_batches", norm_sd=0.5, seed=43):
    """
    Adds a batch effect to the PCA results stored in an AnnData object.

    Batch effects are scaled by the sqrt of sum of component variances to make
    the norm_sd parameter comparable across datasets.

    Parameters:
    - adata: AnnData object containing single-cell data with embeddings in .obsm[layer_embedding]
    - layer_embedding: Name of embedding layer to use (e.g., 'X_pca', 'DM_EigenVectors')
    - batch_col: The column in .obs corresponding to batch information
    - norm_sd: Batch effect strength as fraction of data variance
               (0 = no effect, 1 = effect size equals data std)
    - seed: Random seed for reproducibility

    Returns:
    - Modified AnnData object with added batch effects in .obsm['{layer_embedding}_batch']
    """

    np.random.seed(seed)
    # Extract embeddings
    X_emb = adata.obsm[layer_embedding]

    X_emb_df = pd.DataFrame(X_emb, index=adata.obs_names)

    # Initialize batch-corrected version with the original
    X_emb_batch = X_emb_df.copy()

    # Calculate scaling factor: sqrt of sum of variances across all components
    # This makes batch effect strength comparable across different datasets
    component_vars = np.var(X_emb, axis=0)
    variance_scale = np.sqrt(np.sum(component_vars))

    # Split cells by batch
    cell_batches = adata.obs[batch_col]

    # Add batch effect for each batch
    # norm_sd=0 means no effect, norm_sd=1 means effect ~ data std
    for batch in cell_batches.unique():
        # Generate batch effect scaled by data variance
        batch_effect = np.random.normal(loc=0, scale=norm_sd * variance_scale, size=X_emb.shape[1])
        batch_indices = list(cell_batches[cell_batches == batch].index)
        X_emb_batch.loc[batch_indices] += batch_effect

    # Store the modified embeddings with batch effects
    adata.obsm[f"{layer_embedding}_batch"] = X_emb_batch.values

    return adata


### def add_batch_effect_DM(adata, layer_embedding,batch_col="synth_batches", norm_sd=0.5,seed = 43):