#!/bin/bash

# Check if module command is available and load modules on systems that support it
if command -v module &> /dev/null; then
    module purge
    #module load R/4.3.1-gfbf-2022b
    module load ImageMagick/7.1.0-53-GCCcore-12.2.0 || true
    module load GSL/2.7-GCCcore-12.2.0 || true
    module load cuDNN/8.4.1.50-CUDA-11.7.0 || true
fi

# Set up micromamba environment if needed
eval "$(micromamba shell hook --shell bash 2>/dev/null)" || echo "micromamba not available, assuming environment is already activated"



script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi
cd ${root}/python_method
echo "files are in :$root/python_method" 

#data_dir="$root/data"
#out_dir="$root/benchmark_python"

#echo "data_dir is $data_dir" 

## Run real data ##

data_id=$1
echo "$data_id"
embedding_layer=$2
echo "$embedding_layer"
n_dm=$3
echo "$n_dm"
mode_embedding=$4
echo "$mode_embedding"

# job_number=0
# M2 M3 M4 M5 M6 M7
# 44 45
#$(seq 0.75 0.1 0.85)



if [ "$data_id" == "cluster" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/cluster_data_bm.h5ad
    pop_col="celltype"
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds

    #out_dir = $root/benchmark_python/synthetic/$data_id
# elif [ "$data_id" == "cluster_balanced" ]
#     then
#     data_dir=${root}/data/synthetic/$data_id
#     data_file=${data_dir}/cluster_balanced_anndata.h5ad
#     pops=$(for p in $(seq 1 1 3); do echo M$p; done)
#     #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
#     batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
#     k=30
#     resolution=0.2
#     beta=33
#     downsample=3
#     mem=8g
#     pop_col="celltype"
elif [ "$data_id" == "linear" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/linear_data_bm.h5ad
    pop_col="celltype"
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
elif [ "$data_id" == "branch" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/branch_data_bm.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pop_col="celltype"
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
elif [ "$data_id" == "aging" ]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/aging_hematopoiesis_data_benchmarking.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="midres_celltype_benchmarking"
elif [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/single-cell-atlas-pbmc-sars-cov2_sce.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="cell.type.coarse"
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/bcr_xl_preprocessed_sce.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="cell_type"
elif [[ "$data_id" == "pancreas" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/pancreas_preprocessed_sce.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="Factor.Value.inferred.cell.type...authors.labels."
elif [[ "$data_id" == "levine32" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/levine32_preprocessed_sce.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="cell_type"
fi


#PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC
echo "data_dir is $data_dir" 


if [ -f "$data_file" ]; then
    python data_preprocessing_pipeline.py --file_path "${data_file}" \
        --embedding_layer "${embedding_layer}" \
        --n_dm "${n_dm}" \
        --pop_col "${pop_col}" \
        --mode_embedding "${mode_embedding}" \
        --output_dir "${out_dir}"
    

    python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"

else
     echo "Skipping $data_id - dataset file not found at $data_file"
fi


# if [ "$data_id" == "levine32" ]; then
#     python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
# elif [ "$data_id" == "aging" ]; then
#     python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
# elif [ "$data_id" == "bcr-xl" ]; then
#     python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
# fi

echo "Done"