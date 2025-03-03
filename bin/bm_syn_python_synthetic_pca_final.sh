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
iteration_num=$3
balance_bool=$4
mellon_method=$5
norm_or_not=$6
corrected=$7
hyper=$8
# job_number=0
# M2 M3 M4 M5 M6 M7
# 44 45
#$(seq 0.75 0.1 0.85)



if [ "$data_id" == "cluster" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/synthetic/cluster/cluster_anndata.h5ad
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in mellon meld cna meld_default; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
    pop_col="celltype"
    #out_dir = $root/benchmark_python/synthetic/$data_id
elif [ "$data_id" == "cluster_balanced" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=${data_dir}/cluster_balanced_anndata.h5ad
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in mellon meld cna meld_default; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
    pop_col="celltype"
elif [ "$data_id" == "linear" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/synthetic/linear/linear_anndata.h5ad
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in meld mellon cna meld_default; do echo $m; done)
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
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/synthetic/branch/branch_anndata.h5ad
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    R_methods=$(for m in mellon meld cna meld_default; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=65
    downsample=3
    mem=8g
    pop_col="celltype"
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
                for method in $R_methods
                    do
                    for iteration in $(seq 0 1 $iteration_num)
                        do
                            ((job_number++))
                            if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                                continue
                            fi
                            jobid=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}
                            jobid_2=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}-${mellon_method}-${norm_or_not}-${hyper}-${corrected}

                            save_path=${root}/benchmark_pca/synthetic/$data_id/${jobid_2}
                            save_path_iteration=${root}/benchmark_pca/synthetic/$data_id/${jobid_2}/iteration_${iteration}
                            mkdir -p "$save_path_iteration"
                            echo "Doing $jobid ..."
                            if [[ "$method" == "mellon" ]]; then

                                conda deactivate
                                conda activate DiffAbundance
                                python Mellon_bm_new.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --mode_select centroid \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --n_random_cell 0 \
                                    --mode_embedding PCA \
                                    --n_dm 0 \
                                    --mellon_d_method ${mellon_method} \
                                    --norm_density ${norm_or_not} \
                                    --hyperparameter ${hyper} \
                                    --corrected ${corrected} \
                                    --ls_factor 1.5 \
                                    --ls_mode PCA \
                                    --output_dir $save_path/
                                exit $!
                            
                            elif [[ "$method" == "meld" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python meld_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --mode_select centroid \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --n_random_cell 0 \
                                    --mode_embedding PCA \
                                    --beta $beta \
                                    --k_meld $k \
                                    --output_dir $save_path/
                                exit $!
                            elif [[ "$method" == "meld_default" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python meld_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --mode_select centroid \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --n_random_cell 0 \
                                    --mode_embedding PCA \
                                    --beta 40 \
                                    --k_meld $k \
                                    --output_dir $save_path/
                                exit $!
                            elif [[ "$method" == "cna" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python CNA_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --mode_select centroid \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --n_random_cell 0 \
                                    --mode_embedding PCA \
                                    --k_cna $k \
                                    --output_dir $save_path/
                                exit $!

                            fi                     
                    done					
				done
			done
        done
    done
done

# Submit a slurm array job
jobid="mellon_syn_real_$data_id"
cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
--mem 8g --out '$root/SlurmLog/${jobid_2}_%N_%A_%a.out' --array=1-$job_number \
'$script_path' $1 $2 $3 $4 $5 $6 $7 $8"
echo "$cmd"
eval "$cmd"
