#!/bin/bash

module purge
#module load R/4.3.1-gfbf-2022b
module load ImageMagick/7.1.0-53-GCCcore-12.2.0
module load GSL/2.7-GCCcore-12.2.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
eval "$(conda shell.bash hook)"

# set slurm parameters
time=1-00:00:00
partition=campus-new

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
analysis_layer=$2
balance_bool=$3
mode_embedding=$4
n_dm=$5

# job_number=0
# M2 M3 M4 M5 M6 M7
# 44 45
#$(seq 0.75 0.1 0.85)



if [ "$data_id" == "cluster" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
    pop_col="celltype"
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
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    #0.75 1 1.25 1.5
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=71
    downsample=3
    mem=8g
    pop_col="celltype"
elif [ "$data_id" == "branch" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=65
    downsample=3
    mem=8g
    pop_col="celltype"
elif [ "$data_id" == "aging" ]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for m in CLP Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC; do echo $m; done)
    #Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    #0.75 1 1.25 1.5
    #$(for m in $(seq 0 0.1 1) $(seq 1 0.25 1.75) 20; do echo $m; done)
    #batch_vec=$(for m in $(seq 0 0.1 1); do echo $m; done)
    batch_vec=$(for m in 0.0; do echo $m; done)
    k=30
    resolution=1
    beta=80
    # this beta will be changed after the parameter tuning
    downsample=3
    mem=8g
    pop_col="midres_celltype_benchmarking"
elif [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for m in RBC B PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC; do echo $m; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.5
    beta=25
    downsample=3
    mem=32g
    pop_col="cell.type.coarse"
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for m in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC; do echo $m; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=23
    downsample=10
    mem=32g
    pop_col="cell_type"
elif [[ "$data_id" == "pancreas" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    pops=$(for m in delta_cell alpha_cell gamma_cell acinar_cell beta_cell ductal_cell epsilon_cell; do echo $m; done)
    #R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1.2
    beta=80
    downsample=3
    mem=8g
    pop_col="Factor.Value.inferred.cell.type...authors.labels."
elif [[ "$data_id" == "levine32" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=${data_dir}/${data_id}_${mode_embedding}_${n_dm}.h5ad
    # CD4_T_cells CD8_T_cells Pre_B_cells Mature_B_cells Monocytes Basophils
    pops=$(for m in pDCs CD4_T_cells CD8_T_cells Pre_B_cells Mature_B_cells Monocytes Basophils; do echo $m; done)
    #R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=36
    downsample=25
    mem=96g
    pop_col="cell_type"
fi


#PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC
echo "data_dir is $data_dir" 
job_number=0

for p in $pops;
    do
    for seed in 43 44 45
        do
        for enr in $(seq 0.75 0.1 0.95)
        	do
			for batch_sd_num in $batch_vec
                do
                    ((job_number++))
                    if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                        continue
                    fi
					jobid=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}
					echo "Doing $jobid ..."
                    conda deactivate
                    conda activate DiffAbundance
                    DIRECTORY=$data_dir/${jobid}/
                    mkdir "$DIRECTORY"
                    python generate_bm_data.py \
                        --file_path ${data_file} \
                        --pop ${p} \
                        --pop_enr $enr \
                        --pop_column ${pop_col} \
                        --ds_type $data_id \
                        --batch_sd ${batch_sd_num} \
                        --n_conditions 2 \
                        --n_replicates 3 \
                        --n_batches 2 \
                        --seed ${seed} \
                        --condition_balance 1 \
                        --m 2 \
                        --a_logit 0.5 \
                        --mode_embedding PCA \
                        --layer_embedding X_pca \
                        --balance $balance_bool \
                        --output_dir $DIRECTORY/
                    exit $!					

			done
        done
    done
done

# Submit a slurm array job
jobid="mellon_syn_real_$data_id"
cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
--mem 8g --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \
'$script_path' $1 $2 $3 $4 $5"
echo "$cmd"
eval "$cmd"