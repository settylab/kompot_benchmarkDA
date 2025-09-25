#!/bin/bash

# Note: This script is designed to be user-agnostic and does not use the module system

# Set up micromamba environment if needed
eval "$(micromamba shell hook --shell bash 2>/dev/null)" || echo "micromamba not available, assuming environment is already activated"

# Activate benchmarkda environment
if command -v micromamba &> /dev/null; then
    micromamba activate benchmarkda 2>/dev/null || echo "Failed to activate benchmarkda environment"
elif [ -n "${MAMBA_EXE}" ]; then
    # Use MAMBA_EXE if available (more reliable)
    MAMBA_PREFIX=${MAMBA_EXE%/bin/mamba}
    export PATH="${MAMBA_PREFIX}/envs/benchmarkda/bin:$PATH"
    export CONDA_DEFAULT_ENV=benchmarkda
fi



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
    data_file=${data_dir}/cluster.h5ad
    pop_col="celltype"
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
elif [ "$data_id" == "linear" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/linear.h5ad
    pop_col="celltype"
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
elif [ "$data_id" == "branch" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/branch.h5ad
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
    data_file=${data_dir}/covid19-pbmc.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="cell.type.coarse"
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/bcr-xl.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="cell_type"
elif [[ "$data_id" == "pancreas" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/pancreas.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="Factor.Value.inferred.cell.type...authors.labels."
elif [[ "$data_id" == "levine32" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/levine32.h5ad
    out_dir=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    out_dir_rds=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.rds
    pop_col="cell_type"
fi


#PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC
echo "data_dir is $data_dir" 

# Use MAMBA_EXE to ensure we run in the correct environment
if [ -n "${MAMBA_EXE}" ]; then
    ${MAMBA_EXE} run -n benchmarkda python data_preprocessing_pipeline.py --file_path "${data_file}" \
        --embedding_layer "${embedding_layer}" \
        --n_dm "${n_dm}" \
        --pop_col "${pop_col}" \
        --mode_embedding "${mode_embedding}" \
        --output_dir "${out_dir}"

    ${MAMBA_EXE} run -n benchmarkda python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
else
    python data_preprocessing_pipeline.py --file_path "${data_file}" \
        --embedding_layer "${embedding_layer}" \
        --n_dm "${n_dm}" \
        --pop_col "${pop_col}" \
        --mode_embedding "${mode_embedding}" \
        --output_dir "${out_dir}"

    python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
fi


# if [ "$data_id" == "levine32" ]; then
#     python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
# elif [ "$data_id" == "aging" ]; then
#     python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
# elif [ "$data_id" == "bcr-xl" ]; then
#     python anndata_rds_transfer.py --input_file_path "${out_dir}" --output_file_path "${out_dir_rds}"
# fi

echo "Done"